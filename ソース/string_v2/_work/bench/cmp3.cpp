#include <bits/stdc++.h>
#include <atcoder/string>
using namespace std;
using ll = long long;
template<class T> int sizei(const T& x){return (int)x.size();}
#define rep(i,n) for(ll i=0;i<(ll)(n);++i)
#ifndef GG_STR_BASE
#define GG_STR_BASE
// ---- 文字列基盤（文字列第二世代）  文字の種類 Alphabet・文字の数と位置・ランレングス・変換・比較・区切り ----
// 列（string / vector<整数>）を受け取る関数はどちらでもよい。文字列リテラル（"abc"）も渡せる
// 範囲は半開区間 [l, r)（0 始まり）。位置が「ない」ときは -1。使える範囲の外の引数はメッセージを出して止める
namespace str_internal {
[[noreturn]] inline void fail(const char* func, const std::string& msg) {
    std::cerr << func << ": " << msg << "\n";
    std::abort();
}
// 文字列リテラルを std::string に（それ以外はそのまま参照を返す）
template <class T>
const T& seq(const T& s) { return s; }
template <std::size_t N>
std::string seq(const char (&s)[N]) { return std::string(s); }
inline std::string seq(const char* s) { return std::string(s); }
template <class T>
using seq_t = std::decay_t<decltype(seq(std::declval<const T&>()))>;
template <class T>
using elem_t = std::decay_t<decltype(std::declval<const seq_t<T>&>()[0])>;
// 要素を「値」に（char は unsigned char として 0〜255。負の添字にしない）
template <class T>
long long code_of(T x) {
    if constexpr (std::is_same<T, char>::value || std::is_same<T, signed char>::value || std::is_same<T, unsigned char>::value)
        return (long long)(unsigned char)x;
    else
        return (long long)x;
}
// 範囲の確認（止めるときの文字列作りは別の関数に分けて、確認そのものを軽くする）
[[noreturn]] __attribute__((noinline, cold)) inline void range_fail(const char* func, long long l, long long r, long long n) {
    fail(func, "範囲 [l, r) は 0 ≤ l ≤ r ≤ 長さ で指定する（l = " + std::to_string(l) + ", r = " + std::to_string(r) + ", 長さ = " + std::to_string(n) + "）");
}
inline void check_range(const char* func, long long l, long long r, long long n) {
    if (__builtin_expect(l < 0 || r > n || l > r, 0)) range_fail(func, l, r, n);
}
// a[lo, hi) の中で x 以上の最初の位置（a は昇順）。分岐の少ない二分探索
inline int lower_index(const int* a, int lo, int hi, int x) {
    int len = hi - lo;
    if (len <= 0) return lo;
    const int* b = a + lo;
    while (len > 1) {
        const int half = len >> 1;
        b += (b[half - 1] < x) ? half : 0;
        len -= half;
    }
    return (int)(b - a) + (*b < x);
}
// 昇順の a の中で x より大きい最初の位置（std::upper_bound と同じ。デバッグビルドでも O(log N)）
template <class T>
int upper_index(const std::vector<T>& a, const T& x) {
    int lo = 0, hi = (int)a.size();
    while (lo < hi) {
        const int mid = (int)((unsigned)(lo + hi) >> 1);
        if (!(x < a[mid])) lo = mid + 1;
        else hi = mid;
    }
    return lo;
}
// 数え上げの足し算：T が整数型なら溢れたら止める（mint・double などはそのまま）
template <class T>
T checked_add(const T& a, const T& b, const char* func) {
    if constexpr (std::is_integral<T>::value) {
        T r;
        if (__builtin_add_overflow(a, b, &r)) fail(func, "値が型に収まらない（答えを mod で求めるなら math_v2 の ModInt998 などを使う）");
        return r;
    } else {
        return a + b;
    }
}
// 数え上げの途中の足し算：T が整数型なら溢れたことを印 oa に残すだけ（途中で溢れても答えに使わなければ止めない。最後に印を見て止める）
template <class T>
void add_count(T& a, bool& oa, const T& b, bool ob) {
    if constexpr (std::is_integral<T>::value) {
        T r;
        oa = __builtin_add_overflow(a, b, &r) || oa || ob;
        a = r;
    } else {
        a = a + b;
    }
}
[[noreturn]] inline void count_overflow(const char* func) { fail(func, "値が型に収まらない（答えを mod で求めるなら math_v2 の ModInt998 などを使う）"); }
// 列 → 値の列（char は 0〜255）
template <class Seq>
std::vector<long long> codes_of(const Seq& s0) {
    const auto& s = seq(s0);
    std::vector<long long> c(s.size());
    for (std::size_t i = 0; i < s.size(); ++i) c[i] = code_of(s[i]);
    return c;
}
// 辞書順の比較（-1 / 0 / 1）。char は unsigned char として比べる（std::string と同じ順）
template <class A, class B>
int lex_compare(const A& a, const B& b) {
    const std::size_t n = a.size(), m = b.size(), k = n < m ? n : m;
    for (std::size_t i = 0; i < k; ++i) {
        const long long x = code_of(a[i]), y = code_of(b[i]);
        if (x != y) return x < y ? -1 : 1;
    }
    return n == m ? 0 : (n < m ? -1 : 1);
}
}  // namespace str_internal

// ---- 部分文字列の位置 [l, r) ----
struct StrSpan {
    int l = 0, r = 0;
    int len() const { return r - l; }
    bool empty() const { return r <= l; }
    friend bool operator==(const StrSpan& a, const StrSpan& b) { return a.l == b.l && a.r == b.r; }
    friend bool operator!=(const StrSpan& a, const StrSpan& b) { return !(a == b); }
    friend std::ostream& operator<<(std::ostream& os, const StrSpan& sp) { return os << "[" << sp.l << ", " << sp.r << ")"; }  // デバッグ用
};
template <class Seq>
str_internal::seq_t<Seq> substr_of(const Seq& s0, StrSpan sp) {  // s[sp.l, sp.r) のコピー（string なら string、vector なら vector）
    const auto& s = str_internal::seq(s0);
    str_internal::check_range("substr_of", sp.l, sp.r, (long long)s.size());
    return str_internal::seq_t<Seq>(s.begin() + sp.l, s.begin() + sp.r);
}

// ---- 文字の種類（Alphabet）：文字・値 → 番号 0〜size()-1。番号の順が辞書順として扱われる ----
//   Alphabet::lower()（a〜z。既定）/ upper()（A〜Z）/ digits()（0〜9）/ letters()（A〜Z, a〜z の 52 種。ASCII 順）
//   alnum()（0〜9, A〜Z, a〜z の 62 種。ASCII 順）/ ascii()（0〜127）/ bytes()（0〜255。日本語など 128 以上の char も）
//   of("()[]")（並べた文字だけ。並べた順が番号の順 = 辞書順）/ ints(K)（整数 0〜K-1）/ ints(lo, hi)（整数 [lo, hi)）
//   範囲外の文字・値を id に渡すとメッセージを出して止まる（contains で先に確かめられる）
struct Alphabet {
    Alphabet() : Alphabet(of("abcdefghijklmnopqrstuvwxyz")) {}
    static Alphabet lower() { return of("abcdefghijklmnopqrstuvwxyz"); }
    static Alphabet upper() { return of("ABCDEFGHIJKLMNOPQRSTUVWXYZ"); }
    static Alphabet digits() { return of("0123456789"); }
    static Alphabet letters() { return of("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"); }
    static Alphabet alnum() { return of("0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"); }
    static Alphabet ascii() { return code_range(0, 128); }
    static Alphabet bytes() { return code_range(0, 256); }
    static Alphabet of(const std::string& cs) {  // 並べた文字（同じ文字を 2 回並べると止まる）
        if (cs.empty()) str_internal::fail("Alphabet::of", "文字が 1 つもない");
        Alphabet A(0);
        A.kind_ = CHARS;
        for (std::size_t i = 0; i < cs.size(); ++i) {
            const unsigned char c = (unsigned char)cs[i];
            if (A.table_[c] >= 0) str_internal::fail("Alphabet::of", std::string("同じ文字が 2 回ある: '") + cs[i] + "'");
            A.table_[c] = (short)i;
            A.chars_[i] = c;
        }
        A.size_ = (int)cs.size();
        return A;
    }
    static Alphabet ints(long long K) { return ints(0, K); }  // 整数 0〜K-1
    static Alphabet ints(long long lo, long long hi) {         // 整数 [lo, hi)
        if (hi <= lo || (unsigned long long)hi - (unsigned long long)lo > 1000000000ULL)
            str_internal::fail("Alphabet::ints", "lo < hi かつ hi - lo ≤ 10^9 で指定する（大きな値は compress_values で 0〜K-1 に）");
        Alphabet A(0);
        A.kind_ = INTS;
        A.lo_ = lo;
        A.size_ = (int)(hi - lo);
        return A;
    }

    int size() const { return size_; }
    bool is_chars() const { return kind_ == CHARS; }  // 文字の表で決めたものか（of・lower など）。false なら ints
    template <class T>
    bool contains(T x) const { return index_of(str_internal::code_of(x)) >= 0; }
    template <class T>
    int id(T x) const {  // 文字・値の番号（範囲外は止める）
        const long long v = str_internal::code_of(x);
        const int r = index_of(v);
        if (r < 0) out_of_range(v);
        return r;
    }
    long long value(int i) const {  // 番号 i の値（文字なら文字コード 0〜255）
        check_index(i);
        return kind_ == CHARS ? (long long)chars_[i] : lo_ + i;
    }
    char ch(int i) const {  // 番号 i の文字（ints で 0〜255 の外の値なら止める）
        const long long v = value(i);
        if (v < 0 || v > 255) str_internal::fail("Alphabet::ch", "番号 " + std::to_string(i) + " の値 " + std::to_string(v) + " は文字にできない");
        return (char)(unsigned char)v;
    }
    std::string chars() const {  // 文字を番号の順に並べたもの（ints なら 0〜255 の範囲のときだけ）
        std::string s;
        for (int i = 0; i < size_; ++i) s += ch(i);
        return s;
    }

   private:
    enum Kind : unsigned char { CHARS, INTS };
    explicit Alphabet(int) { table_.fill(-1); }
    static Alphabet code_range(int lo, int hi) {
        Alphabet A(0);
        A.kind_ = CHARS;
        for (int c = lo; c < hi; ++c) A.table_[c] = (short)(c - lo), A.chars_[c - lo] = (unsigned char)c;
        A.size_ = hi - lo;
        return A;
    }
    int index_of(long long v) const {
        if (kind_ == CHARS) return (unsigned long long)v < 256ULL ? table_[(std::size_t)v] : -1;
        const unsigned long long d = (unsigned long long)v - (unsigned long long)lo_;
        return d < (unsigned long long)size_ ? (int)d : -1;
    }
    void check_index(int i) const {
        if (i < 0 || i >= size_) str_internal::fail("Alphabet", "番号 " + std::to_string(i) + " は 0〜" + std::to_string(size_ - 1) + " の外");
    }
    [[noreturn]] void out_of_range(long long v) const {
        std::string what = "値 " + std::to_string(v) + " は";
        if (kind_ == CHARS && v >= 32 && v < 127) what = std::string("文字 '") + (char)v + "'（コード " + std::to_string(v) + "）は";
        str_internal::fail("Alphabet::id", what + "文字の種類に入っていない（種類の数 " + std::to_string(size_) + "）");
    }
    Kind kind_ = CHARS;
    int size_ = 0;
    long long lo_ = 0;
    std::array<short, 256> table_{};
    std::array<unsigned char, 256> chars_{};
};

template <class Seq>
std::vector<int> to_ids(const Seq& s0, const Alphabet& A = Alphabet::lower()) {  // 各文字の番号の列
    const auto& s = str_internal::seq(s0);
    std::vector<int> r(s.size());
    for (std::size_t i = 0; i < s.size(); ++i) r[i] = A.id(s[i]);
    return r;
}
inline std::string from_ids(const std::vector<int>& ids, const Alphabet& A = Alphabet::lower()) {  // 番号の列 → 文字列
    std::string s(ids.size(), '\0');
    for (std::size_t i = 0; i < ids.size(); ++i) s[i] = A.ch(ids[i]);
    return s;
}

// ---- 数の列を 0〜K-1 に圧縮（大小の順を保つ）：Alphabet::ints(K) と組み合わせる ----
struct CompressedValues {
    std::vector<int> ids;            // 各要素の番号（0〜K-1。値の小さい順）
    std::vector<long long> values;   // 番号 → 元の値（昇順。K = values.size()）
    int size() const { return (int)values.size(); }
    Alphabet alphabet() const { return Alphabet::ints(values.empty() ? 1 : (long long)values.size()); }
};
template <class Seq>
CompressedValues compress_values(const Seq& s0) {
    const auto& s = str_internal::seq(s0);
    const std::size_t n = s.size();
    std::vector<std::pair<long long, int>> vi(n);
    for (std::size_t i = 0; i < n; ++i) vi[i] = {str_internal::code_of(s[i]), (int)i};
    std::sort(vi.begin(), vi.end());
    CompressedValues c;
    c.ids.resize(n);
    for (std::size_t i = 0; i < n; ++i) {
        if (i == 0 || vi[i].first != vi[i - 1].first) c.values.push_back(vi[i].first);
        c.ids[vi[i].second] = (int)c.values.size() - 1;
    }
    return c;
}

// ---- 文字の数 ----
template <class Seq>
std::vector<int> count_chars(const Seq& s0, const Alphabet& A = Alphabet::lower()) {  // 番号ごとの個数（長さ A.size()）
    const auto& s = str_internal::seq(s0);
    std::vector<int> c(A.size(), 0);
    for (std::size_t i = 0; i < s.size(); ++i) ++c[A.id(s[i])];
    return c;
}

// ---- 区間の文字の数・次 / 前の出現を O(1) で（前処理 O(Nσ)・メモリ (N+1)σ 個の int。N = 2×10^5, σ = 26 で約 21MB） ----
//   メモリを抑えたいとき・σ が大きいときは CharPositions（O(log N)・メモリ O(N + σ)）
struct CharPrefix {
    CharPrefix() = default;
    template <class Seq>
    explicit CharPrefix(const Seq& s0, const Alphabet& A = Alphabet::lower()) : A_(A) {
        const auto& s = str_internal::seq(s0);
        n_ = (int)s.size(), sg_ = A_.size();
        if ((long long)(n_ + 1) * sg_ > 600000000LL) str_internal::fail("CharPrefix", "(N + 1)σ が 6×10^8 を超える（メモリが足りない。CharPositions を使う）");
        cnt_.assign((std::size_t)(n_ + 1) * sg_, 0);
        std::vector<int> id(n_);
        start_.assign(sg_ + 1, 0);
        for (int i = 0; i < n_; ++i) {
            const int* src = cnt_.data() + (std::size_t)i * sg_;
            int* dst = cnt_.data() + (std::size_t)(i + 1) * sg_;
            std::copy(src, src + sg_, dst);
            id[i] = A_.id(s[i]);
            ++dst[id[i]], ++start_[id[i] + 1];
        }
        for (int k = 0; k < sg_; ++k) start_[k + 1] += start_[k];
        pos_.resize(n_);
        std::vector<int> cur(start_.begin(), start_.end() - 1);
        for (int i = 0; i < n_; ++i) pos_[cur[id[i]]++] = i;
    }
    int size() const { return n_; }
    const Alphabet& alphabet() const { return A_; }
    template <class C>
    int count(int l, int r, C c) const {  // [l, r) にある c の数  O(1)
        str_internal::check_range("CharPrefix::count", l, r, n_);
        const int k = A_.id(c);
        return cnt_[(std::size_t)r * sg_ + k] - cnt_[(std::size_t)l * sg_ + k];
    }
    template <class C>
    bool has(int l, int r, C c) const { return count(l, r, c) > 0; }  // [l, r) に c があるか  O(1)
    std::vector<int> counts(int l, int r) const {  // [l, r) の番号ごとの数  O(σ)
        str_internal::check_range("CharPrefix::counts", l, r, n_);
        std::vector<int> c(sg_);
        for (int k = 0; k < sg_; ++k) c[k] = cnt_[(std::size_t)r * sg_ + k] - cnt_[(std::size_t)l * sg_ + k];
        return c;
    }
    int distinct(int l, int r) const {  // [l, r) にある文字の種類の数  O(σ)
        str_internal::check_range("CharPrefix::distinct", l, r, n_);
        int d = 0;
        for (int k = 0; k < sg_; ++k) d += cnt_[(std::size_t)r * sg_ + k] != cnt_[(std::size_t)l * sg_ + k];
        return d;
    }
    template <class C>
    int next(int i, C c) const {  // i 以上で最初の c の位置（なければ -1。0 ≤ i ≤ N）  O(1)
        check_pos("CharPrefix::next", i);
        const int k = A_.id(c), m = cnt_[(std::size_t)i * sg_ + k];
        return start_[k] + m < start_[k + 1] ? pos_[start_[k] + m] : -1;
    }
    template <class C>
    int prev(int i, C c) const {  // i より前（[0, i)）で最後の c の位置（なければ -1。0 ≤ i ≤ N）  O(1)
        check_pos("CharPrefix::prev", i);
        const int k = A_.id(c), m = cnt_[(std::size_t)i * sg_ + k];
        return m > 0 ? pos_[start_[k] + m - 1] : -1;
    }
    template <class C>
    int kth(C c, long long k) const {  // c の k 番目（0 から）の位置（なければ -1）  O(1)
        const int id = A_.id(c);
        return (k >= 0 && k < start_[id + 1] - start_[id]) ? pos_[start_[id] + k] : -1;
    }
    template <class C>
    int total(C c) const {  // 全体での c の数  O(1)
        const int k = A_.id(c);
        return start_[k + 1] - start_[k];
    }

   private:
    void check_pos(const char* func, int i) const {
        if (i < 0 || i > n_) str_internal::fail(func, "位置 " + std::to_string(i) + " は 0〜N（= " + std::to_string(n_) + "）の外");
    }
    Alphabet A_;
    int n_ = 0, sg_ = 0;
    std::vector<int> cnt_;           // cnt_[i * σ + k] = [0, i) にある番号 k の数
    std::vector<int> start_, pos_;   // 番号ごとの出現位置（pos_[start_[k] .. start_[k+1])）
};

// ---- 文字の出現位置（前処理 O(N + σ)・メモリ O(N + σ)。σ が大きくても使える） ----
struct CharPositions {
    CharPositions() = default;
    template <class Seq>
    explicit CharPositions(const Seq& s0, const Alphabet& A = Alphabet::lower()) : A_(A) {
        const auto& s = str_internal::seq(s0);
        n_ = (int)s.size();
        const int sg = A_.size();
        std::vector<int> id(n_);
        start_.assign(sg + 1, 0);
        for (int i = 0; i < n_; ++i) id[i] = A_.id(s[i]), ++start_[id[i] + 1];
        for (int k = 0; k < sg; ++k) start_[k + 1] += start_[k];
        pos_.resize(n_), occ_.resize(n_);
        std::vector<int> cur(start_.begin(), start_.end() - 1);
        for (int i = 0; i < n_; ++i) occ_[i] = cur[id[i]] - start_[id[i]], pos_[cur[id[i]]++] = i;
    }
    int size() const { return n_; }
    const Alphabet& alphabet() const { return A_; }
    template <class C>
    int next(int i, C c) const {  // i 以上で最初の c の位置（なければ -1。0 ≤ i ≤ N）  O(log N)
        check_pos("CharPositions::next", i);
        const int k = A_.id(c), j = str_internal::lower_index(pos_.data(), start_[k], start_[k + 1], i);
        return j < start_[k + 1] ? pos_[j] : -1;
    }
    template <class C>
    int prev(int i, C c) const {  // i より前（[0, i)）で最後の c の位置（なければ -1。0 ≤ i ≤ N）  O(log N)
        check_pos("CharPositions::prev", i);
        const int k = A_.id(c), j = str_internal::lower_index(pos_.data(), start_[k], start_[k + 1], i);
        return j > start_[k] ? pos_[j - 1] : -1;
    }
    template <class C>
    int kth(C c, long long k) const {  // c の k 番目（0 から）の位置（なければ -1）  O(1)
        const int id = A_.id(c);
        return (k >= 0 && k < start_[id + 1] - start_[id]) ? pos_[start_[id] + k] : -1;
    }
    template <class C>
    int count(int l, int r, C c) const {  // [l, r) にある c の数  O(log N)
        str_internal::check_range("CharPositions::count", l, r, n_);
        const int k = A_.id(c);
        return str_internal::lower_index(pos_.data(), start_[k], start_[k + 1], r) - str_internal::lower_index(pos_.data(), start_[k], start_[k + 1], l);
    }
    template <class C>
    bool has(int l, int r, C c) const { return count(l, r, c) > 0; }
    template <class C>
    int total(C c) const {  // 全体での c の数  O(1)
        const int k = A_.id(c);
        return start_[k + 1] - start_[k];
    }
    template <class C>
    std::vector<int> positions(C c) const {  // c の位置の一覧（昇順）
        const int k = A_.id(c);
        return std::vector<int>(pos_.begin() + start_[k], pos_.begin() + start_[k + 1]);
    }
    int occurrence_index(int i) const {  // s[i] が「その文字の何番目（0 から）の出現」か  O(1)
        if (i < 0 || i >= n_) str_internal::fail("CharPositions::occurrence_index", "位置 " + std::to_string(i) + " は 0〜N-1 の外");
        return occ_[i];
    }

   private:
    void check_pos(const char* func, int i) const {
        if (i < 0 || i > n_) str_internal::fail(func, "位置 " + std::to_string(i) + " は 0〜N（= " + std::to_string(n_) + "）の外");
    }
    Alphabet A_;
    int n_ = 0;
    std::vector<int> start_, pos_, occ_;
};

// ---- ランレングス圧縮：連続する同じ値を {value, len, start} にまとめる ----
template <class T>
struct RunBlock {
    T value{};
    long long len = 0;    // 続く長さ
    long long start = 0;  // 元の列での開始位置
    friend std::ostream& operator<<(std::ostream& os, const RunBlock& b) { return os << "{" << b.value << ", " << b.len << ", " << b.start << "}"; }  // デバッグ用
};
template <class Seq>
std::vector<RunBlock<str_internal::elem_t<Seq>>> run_length(const Seq& s0) {
    const auto& s = str_internal::seq(s0);
    std::vector<RunBlock<str_internal::elem_t<Seq>>> r;
    std::size_t cnt = s.empty() ? 0 : 1;
    for (std::size_t i = 1; i < s.size(); ++i) cnt += !(s[i] == s[i - 1]);
    r.reserve(cnt);
    for (std::size_t i = 0, n = s.size(); i < n;) {
        std::size_t j = i + 1;
        while (j < n && s[j] == s[i]) ++j;
        r.push_back({s[i], (long long)(j - i), (long long)i});
        i = j;
    }
    return r;
}
template <class T>
auto run_length_decode(const std::vector<RunBlock<T>>& runs) {  // 元に戻す（T = char なら string、それ以外は vector<T>）
    using R = std::conditional_t<std::is_same<T, char>::value, std::string, std::vector<T>>;
    long long total = 0;
    for (const auto& b : runs) {
        if (b.len < 0) str_internal::fail("run_length_decode", "len が負");
        total += b.len;
        if (total > 2000000000LL) str_internal::fail("run_length_decode", "元の長さが 2×10^9 を超える");
    }
    R r;
    r.reserve((std::size_t)total);
    for (const auto& b : runs) r.insert(r.end(), (std::size_t)b.len, b.value);
    return r;
}

// ---- 変換 ----
inline std::string to_lower(std::string s) {  // 英大文字を小文字に（ほかはそのまま）
    for (char& c : s)
        if ('A' <= c && c <= 'Z') c = (char)(c - 'A' + 'a');
    return s;
}
inline std::string to_upper(std::string s) {
    for (char& c : s)
        if ('a' <= c && c <= 'z') c = (char)(c - 'a' + 'A');
    return s;
}
inline std::string swap_case(std::string s) {  // 大文字と小文字を入れ替える
    for (char& c : s) {
        if ('A' <= c && c <= 'Z') c = (char)(c - 'A' + 'a');
        else if ('a' <= c && c <= 'z') c = (char)(c - 'a' + 'A');
    }
    return s;
}
template <class Seq>
str_internal::seq_t<Seq> caesar_shift(const Seq& s0, long long k, const Alphabet& A = Alphabet::lower()) {  // A の文字を番号で k ずらす（負も。A にない文字はそのまま）
    auto s = str_internal::seq_t<Seq>(str_internal::seq(s0));
    const long long sg = A.size(), sh = ((k % sg) + sg) % sg;
    for (auto& c : s)
        if (A.contains(c)) c = (std::decay_t<decltype(c)>)A.value((int)((A.id(c) + sh) % sg));
    return s;
}
// 文字の置き換えを重ねる：replace(c, d) で「今 c になっている文字を全部 d に」。apply(s) で s に当てはめる（ABC342 C）
struct CharMap {
    explicit CharMap(const Alphabet& A = Alphabet::lower()) : A_(A), f_(A.size()) {
        for (int i = 0; i < A_.size(); ++i) f_[i] = i;
    }
    template <class C, class D>
    void replace(C c, D d) {  // 今 c になっているものを全部 d に  O(σ)
        const int x = A_.id(c), y = A_.id(d);
        if (x == y) return;
        for (int& v : f_)
            if (v == x) v = y;
    }
    template <class C, class D>
    void set(C c, D d) { f_[A_.id(c)] = A_.id(d); }  // 元の文字 c の行き先を d にする（ほかは変えない）
    template <class C>
    auto operator()(C c) const {  // 元の文字 c の行き先（A にない文字はそのまま）
        return A_.contains(c) ? (C)A_.value(f_[A_.id(c)]) : c;
    }
    template <class Seq>
    str_internal::seq_t<Seq> apply(const Seq& s0) const {  // s の各文字を行き先に（A にない文字はそのまま）
        auto s = str_internal::seq_t<Seq>(str_internal::seq(s0));
        for (auto& c : s) c = (*this)(c);
        return s;
    }
    CharMap compose(const CharMap& g) const {  // 「この置き換えのあとに g」を 1 つにしたもの
        CharMap h(A_);
        for (int i = 0; i < A_.size(); ++i) h.f_[i] = g.f_[f_[i]];
        return h;
    }
    const std::vector<int>& table() const { return f_; }  // 番号 → 行き先の番号

   private:
    Alphabet A_;
    std::vector<int> f_;
};
template <class Seq>
str_internal::seq_t<Seq> reversed(const Seq& s0) {  // 反転したコピー
    const auto& s = str_internal::seq(s0);
    return str_internal::seq_t<Seq>(s.rbegin(), s.rend());
}
template <class Seq>
str_internal::seq_t<Seq> min_with_reverse(const Seq& s0) {  // s と反転のうち辞書順で小さい方（向きを区別しないときの代表。ABC310 C）
    auto s = str_internal::seq_t<Seq>(str_internal::seq(s0));
    auto t = reversed(s);
    return str_internal::lex_compare(t, s) < 0 ? t : s;
}
template <class Seq>
str_internal::seq_t<Seq> sorted_chars(const Seq& s0) {  // 並べ替えたもの（並べ替えて同じかの比較用）
    auto s = str_internal::seq_t<Seq>(str_internal::seq(s0));
    std::sort(s.begin(), s.end(), [](const auto& a, const auto& b) { return str_internal::code_of(a) < str_internal::code_of(b); });
    return s;
}
template <class Seq>
str_internal::seq_t<Seq> repeat_string(const Seq& s0, long long k) {  // s を k 回つなげる
    const auto& s = str_internal::seq(s0);
    if (k < 0) str_internal::fail("repeat_string", "k ≥ 0 で使う");
    if (!s.empty() && k > 2000000000LL / (long long)s.size()) str_internal::fail("repeat_string", "長さが 2×10^9 を超える");
    str_internal::seq_t<Seq> r;
    r.reserve((std::size_t)(k * (long long)s.size()));
    for (long long i = 0; i < k; ++i) r.insert(r.end(), s.begin(), s.end());
    return r;
}

// ---- 比較 ----
template <class A, class B>
int common_prefix_length(const A& a0, const B& b0) {  // 共通の接頭辞の長さ
    const auto& a = str_internal::seq(a0);
    const auto& b = str_internal::seq(b0);
    const std::size_t k = a.size() < b.size() ? a.size() : b.size();
    std::size_t i = 0;
    while (i < k && str_internal::code_of(a[i]) == str_internal::code_of(b[i])) ++i;
    return (int)i;
}
template <class A, class B>
int common_suffix_length(const A& a0, const B& b0) {  // 共通の接尾辞の長さ
    const auto& a = str_internal::seq(a0);
    const auto& b = str_internal::seq(b0);
    const std::size_t k = a.size() < b.size() ? a.size() : b.size();
    std::size_t i = 0;
    while (i < k && str_internal::code_of(a[a.size() - 1 - i]) == str_internal::code_of(b[b.size() - 1 - i])) ++i;
    return (int)i;
}
template <class A, class B>
int hamming_distance(const A& a0, const B& b0) {  // 同じ位置で違う文字の数（長さが同じこと）
    const auto& a = str_internal::seq(a0);
    const auto& b = str_internal::seq(b0);
    if (a.size() != b.size()) str_internal::fail("hamming_distance", "長さが違う");
    int d = 0;
    for (std::size_t i = 0; i < a.size(); ++i) d += str_internal::code_of(a[i]) != str_internal::code_of(b[i]);
    return d;
}
template <class T, class S>
bool is_subsequence(const T& t0, const S& s0) {  // t が s の部分列（飛び飛びでよい）か  O(|s|)。何度も聞くなら 06_subsequence の NextTable
    const auto& t = str_internal::seq(t0);
    const auto& s = str_internal::seq(s0);
    std::size_t j = 0;
    for (std::size_t i = 0; i < s.size() && j < t.size(); ++i)
        if (str_internal::code_of(s[i]) == str_internal::code_of(t[j])) ++j;
    return j == t.size();
}
template <class Seq>
bool is_palindrome(const Seq& s0) {  // 全体が回文か  O(N)。区間を何度も聞くなら 04_palindrome の Palindromes
    const auto& s = str_internal::seq(s0);
    for (std::size_t i = 0, n = s.size(); i < n / 2; ++i)
        if (!(s[i] == s[n - 1 - i])) return false;
    return true;
}

// ---- 区切り・つなげる ----
inline std::vector<std::string> split(const std::string& s, char sep, bool skip_empty = false) {  // sep で区切る（空の部分も残す。skip_empty で除く）
    std::vector<std::string> r;
    std::size_t b = 0;
    for (std::size_t i = 0; i <= s.size(); ++i)
        if (i == s.size() || s[i] == sep) {
            if (!(skip_empty && i == b)) r.push_back(s.substr(b, i - b));
            b = i + 1;
        }
    return r;
}
template <class V>
std::string join(const std::vector<V>& vs, const std::string& sep = " ") {  // 区切りを挟んでつなげる（数も文字列にする。答えの出力に）
    std::string r;
    for (std::size_t i = 0; i < vs.size(); ++i) {
        if (i) r += sep;
        if constexpr (std::is_same<V, std::string>::value) r += vs[i];
        else if constexpr (std::is_same<V, char>::value) r += vs[i];
        else if constexpr (std::is_integral<V>::value) r += std::to_string(vs[i]);
        else {
            std::ostringstream os;
            os << vs[i];
            r += os.str();
        }
    }
    return r;
}
#endif  // GG_STR_BASE

#ifndef GG_STR_BASE
#define GG_STR_BASE
// ---- 文字列基盤（文字列第二世代）  文字の種類 Alphabet・文字の数と位置・ランレングス・変換・比較・区切り ----
// 列（string / vector<整数>）を受け取る関数はどちらでもよい。文字列リテラル（"abc"）も渡せる
// 範囲は半開区間 [l, r)（0 始まり）。位置が「ない」ときは -1。使える範囲の外の引数はメッセージを出して止める
namespace str_internal {
[[noreturn]] inline void fail(const char* func, const std::string& msg) {
    std::cerr << func << ": " << msg << "\n";
    std::abort();
}
// 文字列リテラルを std::string に（それ以外はそのまま参照を返す）
template <class T>
const T& seq(const T& s) { return s; }
template <std::size_t N>
std::string seq(const char (&s)[N]) { return std::string(s); }
inline std::string seq(const char* s) { return std::string(s); }
template <class T>
using seq_t = std::decay_t<decltype(seq(std::declval<const T&>()))>;
template <class T>
using elem_t = std::decay_t<decltype(std::declval<const seq_t<T>&>()[0])>;
// 要素を「値」に（char は unsigned char として 0〜255。負の添字にしない）
template <class T>
long long code_of(T x) {
    if constexpr (std::is_same<T, char>::value || std::is_same<T, signed char>::value || std::is_same<T, unsigned char>::value)
        return (long long)(unsigned char)x;
    else
        return (long long)x;
}
// 範囲の確認（止めるときの文字列作りは別の関数に分けて、確認そのものを軽くする）
[[noreturn]] __attribute__((noinline, cold)) inline void range_fail(const char* func, long long l, long long r, long long n) {
    fail(func, "範囲 [l, r) は 0 ≤ l ≤ r ≤ 長さ で指定する（l = " + std::to_string(l) + ", r = " + std::to_string(r) + ", 長さ = " + std::to_string(n) + "）");
}
inline void check_range(const char* func, long long l, long long r, long long n) {
    if (__builtin_expect(l < 0 || r > n || l > r, 0)) range_fail(func, l, r, n);
}
// a[lo, hi) の中で x 以上の最初の位置（a は昇順）。分岐の少ない二分探索
inline int lower_index(const int* a, int lo, int hi, int x) {
    int len = hi - lo;
    if (len <= 0) return lo;
    const int* b = a + lo;
    while (len > 1) {
        const int half = len >> 1;
        b += (b[half - 1] < x) ? half : 0;
        len -= half;
    }
    return (int)(b - a) + (*b < x);
}
// 昇順の a の中で x より大きい最初の位置（std::upper_bound と同じ。デバッグビルドでも O(log N)）
template <class T>
int upper_index(const std::vector<T>& a, const T& x) {
    int lo = 0, hi = (int)a.size();
    while (lo < hi) {
        const int mid = (int)((unsigned)(lo + hi) >> 1);
        if (!(x < a[mid])) lo = mid + 1;
        else hi = mid;
    }
    return lo;
}
// 数え上げの足し算：T が整数型なら溢れたら止める（mint・double などはそのまま）
template <class T>
T checked_add(const T& a, const T& b, const char* func) {
    if constexpr (std::is_integral<T>::value) {
        T r;
        if (__builtin_add_overflow(a, b, &r)) fail(func, "値が型に収まらない（答えを mod で求めるなら math_v2 の ModInt998 などを使う）");
        return r;
    } else {
        return a + b;
    }
}
// 数え上げの途中の足し算：T が整数型なら溢れたことを印 oa に残すだけ（途中で溢れても答えに使わなければ止めない。最後に印を見て止める）
template <class T>
void add_count(T& a, bool& oa, const T& b, bool ob) {
    if constexpr (std::is_integral<T>::value) {
        T r;
        oa = __builtin_add_overflow(a, b, &r) || oa || ob;
        a = r;
    } else {
        a = a + b;
    }
}
[[noreturn]] inline void count_overflow(const char* func) { fail(func, "値が型に収まらない（答えを mod で求めるなら math_v2 の ModInt998 などを使う）"); }
// 列 → 値の列（char は 0〜255）
template <class Seq>
std::vector<long long> codes_of(const Seq& s0) {
    const auto& s = seq(s0);
    std::vector<long long> c(s.size());
    for (std::size_t i = 0; i < s.size(); ++i) c[i] = code_of(s[i]);
    return c;
}
// 辞書順の比較（-1 / 0 / 1）。char は unsigned char として比べる（std::string と同じ順）
template <class A, class B>
int lex_compare(const A& a, const B& b) {
    const std::size_t n = a.size(), m = b.size(), k = n < m ? n : m;
    for (std::size_t i = 0; i < k; ++i) {
        const long long x = code_of(a[i]), y = code_of(b[i]);
        if (x != y) return x < y ? -1 : 1;
    }
    return n == m ? 0 : (n < m ? -1 : 1);
}
}  // namespace str_internal

// ---- 部分文字列の位置 [l, r) ----
struct StrSpan {
    int l = 0, r = 0;
    int len() const { return r - l; }
    bool empty() const { return r <= l; }
    friend bool operator==(const StrSpan& a, const StrSpan& b) { return a.l == b.l && a.r == b.r; }
    friend bool operator!=(const StrSpan& a, const StrSpan& b) { return !(a == b); }
    friend std::ostream& operator<<(std::ostream& os, const StrSpan& sp) { return os << "[" << sp.l << ", " << sp.r << ")"; }  // デバッグ用
};
template <class Seq>
str_internal::seq_t<Seq> substr_of(const Seq& s0, StrSpan sp) {  // s[sp.l, sp.r) のコピー（string なら string、vector なら vector）
    const auto& s = str_internal::seq(s0);
    str_internal::check_range("substr_of", sp.l, sp.r, (long long)s.size());
    return str_internal::seq_t<Seq>(s.begin() + sp.l, s.begin() + sp.r);
}

// ---- 文字の種類（Alphabet）：文字・値 → 番号 0〜size()-1。番号の順が辞書順として扱われる ----
//   Alphabet::lower()（a〜z。既定）/ upper()（A〜Z）/ digits()（0〜9）/ letters()（A〜Z, a〜z の 52 種。ASCII 順）
//   alnum()（0〜9, A〜Z, a〜z の 62 種。ASCII 順）/ ascii()（0〜127）/ bytes()（0〜255。日本語など 128 以上の char も）
//   of("()[]")（並べた文字だけ。並べた順が番号の順 = 辞書順）/ ints(K)（整数 0〜K-1）/ ints(lo, hi)（整数 [lo, hi)）
//   範囲外の文字・値を id に渡すとメッセージを出して止まる（contains で先に確かめられる）
struct Alphabet {
    Alphabet() : Alphabet(of("abcdefghijklmnopqrstuvwxyz")) {}
    static Alphabet lower() { return of("abcdefghijklmnopqrstuvwxyz"); }
    static Alphabet upper() { return of("ABCDEFGHIJKLMNOPQRSTUVWXYZ"); }
    static Alphabet digits() { return of("0123456789"); }
    static Alphabet letters() { return of("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"); }
    static Alphabet alnum() { return of("0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"); }
    static Alphabet ascii() { return code_range(0, 128); }
    static Alphabet bytes() { return code_range(0, 256); }
    static Alphabet of(const std::string& cs) {  // 並べた文字（同じ文字を 2 回並べると止まる）
        if (cs.empty()) str_internal::fail("Alphabet::of", "文字が 1 つもない");
        Alphabet A(0);
        A.kind_ = CHARS;
        for (std::size_t i = 0; i < cs.size(); ++i) {
            const unsigned char c = (unsigned char)cs[i];
            if (A.table_[c] >= 0) str_internal::fail("Alphabet::of", std::string("同じ文字が 2 回ある: '") + cs[i] + "'");
            A.table_[c] = (short)i;
            A.chars_[i] = c;
        }
        A.size_ = (int)cs.size();
        return A;
    }
    static Alphabet ints(long long K) { return ints(0, K); }  // 整数 0〜K-1
    static Alphabet ints(long long lo, long long hi) {         // 整数 [lo, hi)
        if (hi <= lo || (unsigned long long)hi - (unsigned long long)lo > 1000000000ULL)
            str_internal::fail("Alphabet::ints", "lo < hi かつ hi - lo ≤ 10^9 で指定する（大きな値は compress_values で 0〜K-1 に）");
        Alphabet A(0);
        A.kind_ = INTS;
        A.lo_ = lo;
        A.size_ = (int)(hi - lo);
        return A;
    }

    int size() const { return size_; }
    bool is_chars() const { return kind_ == CHARS; }  // 文字の表で決めたものか（of・lower など）。false なら ints
    template <class T>
    bool contains(T x) const { return index_of(str_internal::code_of(x)) >= 0; }
    template <class T>
    int id(T x) const {  // 文字・値の番号（範囲外は止める）
        const long long v = str_internal::code_of(x);
        const int r = index_of(v);
        if (r < 0) out_of_range(v);
        return r;
    }
    long long value(int i) const {  // 番号 i の値（文字なら文字コード 0〜255）
        check_index(i);
        return kind_ == CHARS ? (long long)chars_[i] : lo_ + i;
    }
    char ch(int i) const {  // 番号 i の文字（ints で 0〜255 の外の値なら止める）
        const long long v = value(i);
        if (v < 0 || v > 255) str_internal::fail("Alphabet::ch", "番号 " + std::to_string(i) + " の値 " + std::to_string(v) + " は文字にできない");
        return (char)(unsigned char)v;
    }
    std::string chars() const {  // 文字を番号の順に並べたもの（ints なら 0〜255 の範囲のときだけ）
        std::string s;
        for (int i = 0; i < size_; ++i) s += ch(i);
        return s;
    }

   private:
    enum Kind : unsigned char { CHARS, INTS };
    explicit Alphabet(int) { table_.fill(-1); }
    static Alphabet code_range(int lo, int hi) {
        Alphabet A(0);
        A.kind_ = CHARS;
        for (int c = lo; c < hi; ++c) A.table_[c] = (short)(c - lo), A.chars_[c - lo] = (unsigned char)c;
        A.size_ = hi - lo;
        return A;
    }
    int index_of(long long v) const {
        if (kind_ == CHARS) return (unsigned long long)v < 256ULL ? table_[(std::size_t)v] : -1;
        const unsigned long long d = (unsigned long long)v - (unsigned long long)lo_;
        return d < (unsigned long long)size_ ? (int)d : -1;
    }
    void check_index(int i) const {
        if (i < 0 || i >= size_) str_internal::fail("Alphabet", "番号 " + std::to_string(i) + " は 0〜" + std::to_string(size_ - 1) + " の外");
    }
    [[noreturn]] void out_of_range(long long v) const {
        std::string what = "値 " + std::to_string(v) + " は";
        if (kind_ == CHARS && v >= 32 && v < 127) what = std::string("文字 '") + (char)v + "'（コード " + std::to_string(v) + "）は";
        str_internal::fail("Alphabet::id", what + "文字の種類に入っていない（種類の数 " + std::to_string(size_) + "）");
    }
    Kind kind_ = CHARS;
    int size_ = 0;
    long long lo_ = 0;
    std::array<short, 256> table_{};
    std::array<unsigned char, 256> chars_{};
};

template <class Seq>
std::vector<int> to_ids(const Seq& s0, const Alphabet& A = Alphabet::lower()) {  // 各文字の番号の列
    const auto& s = str_internal::seq(s0);
    std::vector<int> r(s.size());
    for (std::size_t i = 0; i < s.size(); ++i) r[i] = A.id(s[i]);
    return r;
}
inline std::string from_ids(const std::vector<int>& ids, const Alphabet& A = Alphabet::lower()) {  // 番号の列 → 文字列
    std::string s(ids.size(), '\0');
    for (std::size_t i = 0; i < ids.size(); ++i) s[i] = A.ch(ids[i]);
    return s;
}

// ---- 数の列を 0〜K-1 に圧縮（大小の順を保つ）：Alphabet::ints(K) と組み合わせる ----
struct CompressedValues {
    std::vector<int> ids;            // 各要素の番号（0〜K-1。値の小さい順）
    std::vector<long long> values;   // 番号 → 元の値（昇順。K = values.size()）
    int size() const { return (int)values.size(); }
    Alphabet alphabet() const { return Alphabet::ints(values.empty() ? 1 : (long long)values.size()); }
};
template <class Seq>
CompressedValues compress_values(const Seq& s0) {
    const auto& s = str_internal::seq(s0);
    const std::size_t n = s.size();
    std::vector<std::pair<long long, int>> vi(n);
    for (std::size_t i = 0; i < n; ++i) vi[i] = {str_internal::code_of(s[i]), (int)i};
    std::sort(vi.begin(), vi.end());
    CompressedValues c;
    c.ids.resize(n);
    for (std::size_t i = 0; i < n; ++i) {
        if (i == 0 || vi[i].first != vi[i - 1].first) c.values.push_back(vi[i].first);
        c.ids[vi[i].second] = (int)c.values.size() - 1;
    }
    return c;
}

// ---- 文字の数 ----
template <class Seq>
std::vector<int> count_chars(const Seq& s0, const Alphabet& A = Alphabet::lower()) {  // 番号ごとの個数（長さ A.size()）
    const auto& s = str_internal::seq(s0);
    std::vector<int> c(A.size(), 0);
    for (std::size_t i = 0; i < s.size(); ++i) ++c[A.id(s[i])];
    return c;
}

// ---- 区間の文字の数・次 / 前の出現を O(1) で（前処理 O(Nσ)・メモリ (N+1)σ 個の int。N = 2×10^5, σ = 26 で約 21MB） ----
//   メモリを抑えたいとき・σ が大きいときは CharPositions（O(log N)・メモリ O(N + σ)）
struct CharPrefix {
    CharPrefix() = default;
    template <class Seq>
    explicit CharPrefix(const Seq& s0, const Alphabet& A = Alphabet::lower()) : A_(A) {
        const auto& s = str_internal::seq(s0);
        n_ = (int)s.size(), sg_ = A_.size();
        if ((long long)(n_ + 1) * sg_ > 600000000LL) str_internal::fail("CharPrefix", "(N + 1)σ が 6×10^8 を超える（メモリが足りない。CharPositions を使う）");
        cnt_.assign((std::size_t)(n_ + 1) * sg_, 0);
        std::vector<int> id(n_);
        start_.assign(sg_ + 1, 0);
        for (int i = 0; i < n_; ++i) {
            const int* src = cnt_.data() + (std::size_t)i * sg_;
            int* dst = cnt_.data() + (std::size_t)(i + 1) * sg_;
            std::copy(src, src + sg_, dst);
            id[i] = A_.id(s[i]);
            ++dst[id[i]], ++start_[id[i] + 1];
        }
        for (int k = 0; k < sg_; ++k) start_[k + 1] += start_[k];
        pos_.resize(n_);
        std::vector<int> cur(start_.begin(), start_.end() - 1);
        for (int i = 0; i < n_; ++i) pos_[cur[id[i]]++] = i;
    }
    int size() const { return n_; }
    const Alphabet& alphabet() const { return A_; }
    template <class C>
    int count(int l, int r, C c) const {  // [l, r) にある c の数  O(1)
        str_internal::check_range("CharPrefix::count", l, r, n_);
        const int k = A_.id(c);
        return cnt_[(std::size_t)r * sg_ + k] - cnt_[(std::size_t)l * sg_ + k];
    }
    template <class C>
    bool has(int l, int r, C c) const { return count(l, r, c) > 0; }  // [l, r) に c があるか  O(1)
    std::vector<int> counts(int l, int r) const {  // [l, r) の番号ごとの数  O(σ)
        str_internal::check_range("CharPrefix::counts", l, r, n_);
        std::vector<int> c(sg_);
        for (int k = 0; k < sg_; ++k) c[k] = cnt_[(std::size_t)r * sg_ + k] - cnt_[(std::size_t)l * sg_ + k];
        return c;
    }
    int distinct(int l, int r) const {  // [l, r) にある文字の種類の数  O(σ)
        str_internal::check_range("CharPrefix::distinct", l, r, n_);
        int d = 0;
        for (int k = 0; k < sg_; ++k) d += cnt_[(std::size_t)r * sg_ + k] != cnt_[(std::size_t)l * sg_ + k];
        return d;
    }
    template <class C>
    int next(int i, C c) const {  // i 以上で最初の c の位置（なければ -1。0 ≤ i ≤ N）  O(1)
        check_pos("CharPrefix::next", i);
        const int k = A_.id(c), m = cnt_[(std::size_t)i * sg_ + k];
        return start_[k] + m < start_[k + 1] ? pos_[start_[k] + m] : -1;
    }
    template <class C>
    int prev(int i, C c) const {  // i より前（[0, i)）で最後の c の位置（なければ -1。0 ≤ i ≤ N）  O(1)
        check_pos("CharPrefix::prev", i);
        const int k = A_.id(c), m = cnt_[(std::size_t)i * sg_ + k];
        return m > 0 ? pos_[start_[k] + m - 1] : -1;
    }
    template <class C>
    int kth(C c, long long k) const {  // c の k 番目（0 から）の位置（なければ -1）  O(1)
        const int id = A_.id(c);
        return (k >= 0 && k < start_[id + 1] - start_[id]) ? pos_[start_[id] + k] : -1;
    }
    template <class C>
    int total(C c) const {  // 全体での c の数  O(1)
        const int k = A_.id(c);
        return start_[k + 1] - start_[k];
    }

   private:
    void check_pos(const char* func, int i) const {
        if (i < 0 || i > n_) str_internal::fail(func, "位置 " + std::to_string(i) + " は 0〜N（= " + std::to_string(n_) + "）の外");
    }
    Alphabet A_;
    int n_ = 0, sg_ = 0;
    std::vector<int> cnt_;           // cnt_[i * σ + k] = [0, i) にある番号 k の数
    std::vector<int> start_, pos_;   // 番号ごとの出現位置（pos_[start_[k] .. start_[k+1])）
};

// ---- 文字の出現位置（前処理 O(N + σ)・メモリ O(N + σ)。σ が大きくても使える） ----
struct CharPositions {
    CharPositions() = default;
    template <class Seq>
    explicit CharPositions(const Seq& s0, const Alphabet& A = Alphabet::lower()) : A_(A) {
        const auto& s = str_internal::seq(s0);
        n_ = (int)s.size();
        const int sg = A_.size();
        std::vector<int> id(n_);
        start_.assign(sg + 1, 0);
        for (int i = 0; i < n_; ++i) id[i] = A_.id(s[i]), ++start_[id[i] + 1];
        for (int k = 0; k < sg; ++k) start_[k + 1] += start_[k];
        pos_.resize(n_), occ_.resize(n_);
        std::vector<int> cur(start_.begin(), start_.end() - 1);
        for (int i = 0; i < n_; ++i) occ_[i] = cur[id[i]] - start_[id[i]], pos_[cur[id[i]]++] = i;
    }
    int size() const { return n_; }
    const Alphabet& alphabet() const { return A_; }
    template <class C>
    int next(int i, C c) const {  // i 以上で最初の c の位置（なければ -1。0 ≤ i ≤ N）  O(log N)
        check_pos("CharPositions::next", i);
        const int k = A_.id(c), j = str_internal::lower_index(pos_.data(), start_[k], start_[k + 1], i);
        return j < start_[k + 1] ? pos_[j] : -1;
    }
    template <class C>
    int prev(int i, C c) const {  // i より前（[0, i)）で最後の c の位置（なければ -1。0 ≤ i ≤ N）  O(log N)
        check_pos("CharPositions::prev", i);
        const int k = A_.id(c), j = str_internal::lower_index(pos_.data(), start_[k], start_[k + 1], i);
        return j > start_[k] ? pos_[j - 1] : -1;
    }
    template <class C>
    int kth(C c, long long k) const {  // c の k 番目（0 から）の位置（なければ -1）  O(1)
        const int id = A_.id(c);
        return (k >= 0 && k < start_[id + 1] - start_[id]) ? pos_[start_[id] + k] : -1;
    }
    template <class C>
    int count(int l, int r, C c) const {  // [l, r) にある c の数  O(log N)
        str_internal::check_range("CharPositions::count", l, r, n_);
        const int k = A_.id(c);
        return str_internal::lower_index(pos_.data(), start_[k], start_[k + 1], r) - str_internal::lower_index(pos_.data(), start_[k], start_[k + 1], l);
    }
    template <class C>
    bool has(int l, int r, C c) const { return count(l, r, c) > 0; }
    template <class C>
    int total(C c) const {  // 全体での c の数  O(1)
        const int k = A_.id(c);
        return start_[k + 1] - start_[k];
    }
    template <class C>
    std::vector<int> positions(C c) const {  // c の位置の一覧（昇順）
        const int k = A_.id(c);
        return std::vector<int>(pos_.begin() + start_[k], pos_.begin() + start_[k + 1]);
    }
    int occurrence_index(int i) const {  // s[i] が「その文字の何番目（0 から）の出現」か  O(1)
        if (i < 0 || i >= n_) str_internal::fail("CharPositions::occurrence_index", "位置 " + std::to_string(i) + " は 0〜N-1 の外");
        return occ_[i];
    }

   private:
    void check_pos(const char* func, int i) const {
        if (i < 0 || i > n_) str_internal::fail(func, "位置 " + std::to_string(i) + " は 0〜N（= " + std::to_string(n_) + "）の外");
    }
    Alphabet A_;
    int n_ = 0;
    std::vector<int> start_, pos_, occ_;
};

// ---- ランレングス圧縮：連続する同じ値を {value, len, start} にまとめる ----
template <class T>
struct RunBlock {
    T value{};
    long long len = 0;    // 続く長さ
    long long start = 0;  // 元の列での開始位置
    friend std::ostream& operator<<(std::ostream& os, const RunBlock& b) { return os << "{" << b.value << ", " << b.len << ", " << b.start << "}"; }  // デバッグ用
};
template <class Seq>
std::vector<RunBlock<str_internal::elem_t<Seq>>> run_length(const Seq& s0) {
    const auto& s = str_internal::seq(s0);
    std::vector<RunBlock<str_internal::elem_t<Seq>>> r;
    std::size_t cnt = s.empty() ? 0 : 1;
    for (std::size_t i = 1; i < s.size(); ++i) cnt += !(s[i] == s[i - 1]);
    r.reserve(cnt);
    for (std::size_t i = 0, n = s.size(); i < n;) {
        std::size_t j = i + 1;
        while (j < n && s[j] == s[i]) ++j;
        r.push_back({s[i], (long long)(j - i), (long long)i});
        i = j;
    }
    return r;
}
template <class T>
auto run_length_decode(const std::vector<RunBlock<T>>& runs) {  // 元に戻す（T = char なら string、それ以外は vector<T>）
    using R = std::conditional_t<std::is_same<T, char>::value, std::string, std::vector<T>>;
    long long total = 0;
    for (const auto& b : runs) {
        if (b.len < 0) str_internal::fail("run_length_decode", "len が負");
        total += b.len;
        if (total > 2000000000LL) str_internal::fail("run_length_decode", "元の長さが 2×10^9 を超える");
    }
    R r;
    r.reserve((std::size_t)total);
    for (const auto& b : runs) r.insert(r.end(), (std::size_t)b.len, b.value);
    return r;
}

// ---- 変換 ----
inline std::string to_lower(std::string s) {  // 英大文字を小文字に（ほかはそのまま）
    for (char& c : s)
        if ('A' <= c && c <= 'Z') c = (char)(c - 'A' + 'a');
    return s;
}
inline std::string to_upper(std::string s) {
    for (char& c : s)
        if ('a' <= c && c <= 'z') c = (char)(c - 'a' + 'A');
    return s;
}
inline std::string swap_case(std::string s) {  // 大文字と小文字を入れ替える
    for (char& c : s) {
        if ('A' <= c && c <= 'Z') c = (char)(c - 'A' + 'a');
        else if ('a' <= c && c <= 'z') c = (char)(c - 'a' + 'A');
    }
    return s;
}
template <class Seq>
str_internal::seq_t<Seq> caesar_shift(const Seq& s0, long long k, const Alphabet& A = Alphabet::lower()) {  // A の文字を番号で k ずらす（負も。A にない文字はそのまま）
    auto s = str_internal::seq_t<Seq>(str_internal::seq(s0));
    const long long sg = A.size(), sh = ((k % sg) + sg) % sg;
    for (auto& c : s)
        if (A.contains(c)) c = (std::decay_t<decltype(c)>)A.value((int)((A.id(c) + sh) % sg));
    return s;
}
// 文字の置き換えを重ねる：replace(c, d) で「今 c になっている文字を全部 d に」。apply(s) で s に当てはめる（ABC342 C）
struct CharMap {
    explicit CharMap(const Alphabet& A = Alphabet::lower()) : A_(A), f_(A.size()) {
        for (int i = 0; i < A_.size(); ++i) f_[i] = i;
    }
    template <class C, class D>
    void replace(C c, D d) {  // 今 c になっているものを全部 d に  O(σ)
        const int x = A_.id(c), y = A_.id(d);
        if (x == y) return;
        for (int& v : f_)
            if (v == x) v = y;
    }
    template <class C, class D>
    void set(C c, D d) { f_[A_.id(c)] = A_.id(d); }  // 元の文字 c の行き先を d にする（ほかは変えない）
    template <class C>
    auto operator()(C c) const {  // 元の文字 c の行き先（A にない文字はそのまま）
        return A_.contains(c) ? (C)A_.value(f_[A_.id(c)]) : c;
    }
    template <class Seq>
    str_internal::seq_t<Seq> apply(const Seq& s0) const {  // s の各文字を行き先に（A にない文字はそのまま）
        auto s = str_internal::seq_t<Seq>(str_internal::seq(s0));
        for (auto& c : s) c = (*this)(c);
        return s;
    }
    CharMap compose(const CharMap& g) const {  // 「この置き換えのあとに g」を 1 つにしたもの
        CharMap h(A_);
        for (int i = 0; i < A_.size(); ++i) h.f_[i] = g.f_[f_[i]];
        return h;
    }
    const std::vector<int>& table() const { return f_; }  // 番号 → 行き先の番号

   private:
    Alphabet A_;
    std::vector<int> f_;
};
template <class Seq>
str_internal::seq_t<Seq> reversed(const Seq& s0) {  // 反転したコピー
    const auto& s = str_internal::seq(s0);
    return str_internal::seq_t<Seq>(s.rbegin(), s.rend());
}
template <class Seq>
str_internal::seq_t<Seq> min_with_reverse(const Seq& s0) {  // s と反転のうち辞書順で小さい方（向きを区別しないときの代表。ABC310 C）
    auto s = str_internal::seq_t<Seq>(str_internal::seq(s0));
    auto t = reversed(s);
    return str_internal::lex_compare(t, s) < 0 ? t : s;
}
template <class Seq>
str_internal::seq_t<Seq> sorted_chars(const Seq& s0) {  // 並べ替えたもの（並べ替えて同じかの比較用）
    auto s = str_internal::seq_t<Seq>(str_internal::seq(s0));
    std::sort(s.begin(), s.end(), [](const auto& a, const auto& b) { return str_internal::code_of(a) < str_internal::code_of(b); });
    return s;
}
template <class Seq>
str_internal::seq_t<Seq> repeat_string(const Seq& s0, long long k) {  // s を k 回つなげる
    const auto& s = str_internal::seq(s0);
    if (k < 0) str_internal::fail("repeat_string", "k ≥ 0 で使う");
    if (!s.empty() && k > 2000000000LL / (long long)s.size()) str_internal::fail("repeat_string", "長さが 2×10^9 を超える");
    str_internal::seq_t<Seq> r;
    r.reserve((std::size_t)(k * (long long)s.size()));
    for (long long i = 0; i < k; ++i) r.insert(r.end(), s.begin(), s.end());
    return r;
}

// ---- 比較 ----
template <class A, class B>
int common_prefix_length(const A& a0, const B& b0) {  // 共通の接頭辞の長さ
    const auto& a = str_internal::seq(a0);
    const auto& b = str_internal::seq(b0);
    const std::size_t k = a.size() < b.size() ? a.size() : b.size();
    std::size_t i = 0;
    while (i < k && str_internal::code_of(a[i]) == str_internal::code_of(b[i])) ++i;
    return (int)i;
}
template <class A, class B>
int common_suffix_length(const A& a0, const B& b0) {  // 共通の接尾辞の長さ
    const auto& a = str_internal::seq(a0);
    const auto& b = str_internal::seq(b0);
    const std::size_t k = a.size() < b.size() ? a.size() : b.size();
    std::size_t i = 0;
    while (i < k && str_internal::code_of(a[a.size() - 1 - i]) == str_internal::code_of(b[b.size() - 1 - i])) ++i;
    return (int)i;
}
template <class A, class B>
int hamming_distance(const A& a0, const B& b0) {  // 同じ位置で違う文字の数（長さが同じこと）
    const auto& a = str_internal::seq(a0);
    const auto& b = str_internal::seq(b0);
    if (a.size() != b.size()) str_internal::fail("hamming_distance", "長さが違う");
    int d = 0;
    for (std::size_t i = 0; i < a.size(); ++i) d += str_internal::code_of(a[i]) != str_internal::code_of(b[i]);
    return d;
}
template <class T, class S>
bool is_subsequence(const T& t0, const S& s0) {  // t が s の部分列（飛び飛びでよい）か  O(|s|)。何度も聞くなら 06_subsequence の NextTable
    const auto& t = str_internal::seq(t0);
    const auto& s = str_internal::seq(s0);
    std::size_t j = 0;
    for (std::size_t i = 0; i < s.size() && j < t.size(); ++i)
        if (str_internal::code_of(s[i]) == str_internal::code_of(t[j])) ++j;
    return j == t.size();
}
template <class Seq>
bool is_palindrome(const Seq& s0) {  // 全体が回文か  O(N)。区間を何度も聞くなら 04_palindrome の Palindromes
    const auto& s = str_internal::seq(s0);
    for (std::size_t i = 0, n = s.size(); i < n / 2; ++i)
        if (!(s[i] == s[n - 1 - i])) return false;
    return true;
}

// ---- 区切り・つなげる ----
inline std::vector<std::string> split(const std::string& s, char sep, bool skip_empty = false) {  // sep で区切る（空の部分も残す。skip_empty で除く）
    std::vector<std::string> r;
    std::size_t b = 0;
    for (std::size_t i = 0; i <= s.size(); ++i)
        if (i == s.size() || s[i] == sep) {
            if (!(skip_empty && i == b)) r.push_back(s.substr(b, i - b));
            b = i + 1;
        }
    return r;
}
template <class V>
std::string join(const std::vector<V>& vs, const std::string& sep = " ") {  // 区切りを挟んでつなげる（数も文字列にする。答えの出力に）
    std::string r;
    for (std::size_t i = 0; i < vs.size(); ++i) {
        if (i) r += sep;
        if constexpr (std::is_same<V, std::string>::value) r += vs[i];
        else if constexpr (std::is_same<V, char>::value) r += vs[i];
        else if constexpr (std::is_integral<V>::value) r += std::to_string(vs[i]);
        else {
            std::ostringstream os;
            os << vs[i];
            r += os.str();
        }
    }
    return r;
}
#endif  // GG_STR_BASE
#ifndef GG_ROLLING_HASH
#define GG_ROLLING_HASH
// ---- ローリングハッシュ（文字列第二世代）  mod 2^61 - 1。基数は実行ごとの乱数で、すべての文字列で共通 ----
// 別の文字列の部分文字列どうし・つなげたものも比べられる。累乗の表は基数を変えずに伸ばす
// 1 回の比較で衝突する確率は約 N / 2^61（N は長さ）
namespace str_internal {
__extension__ typedef unsigned __int128 rh_u128;
struct RHCore {
    static constexpr unsigned long long MOD = (1ULL << 61) - 1;
    static unsigned long long mul(unsigned long long a, unsigned long long b) {
        const rh_u128 t = (rh_u128)a * b;
        const unsigned long long r = (unsigned long long)(t >> 61) + ((unsigned long long)t & MOD);
        return r >= MOD ? r - MOD : r;
    }
    static unsigned long long add(unsigned long long a, unsigned long long b) {
        const unsigned long long r = a + b;
        return r >= MOD ? r - MOD : r;
    }
    static unsigned long long sub(unsigned long long a, unsigned long long b) { return a >= b ? a - b : a + MOD - b; }
    static unsigned long long splitmix(unsigned long long x) {
        x += 0x9E3779B97F4A7C15ULL;
        x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
        x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
        return x ^ (x >> 31);
    }
    inline static unsigned long long seed = 0, base1 = 0, base2 = 0;
    inline static std::vector<unsigned long long> pw1, pw2;
    static void init() {
        if (base1) return;
        std::random_device rd;
        seed = splitmix(((unsigned long long)rd() << 32) ^ rd() ^ (unsigned long long)std::chrono::steady_clock::now().time_since_epoch().count());
        const unsigned long long lo = 1ULL << 40, w = MOD - 2 * lo;
        base1 = lo + splitmix(seed ^ 1) % w;
        do base2 = lo + splitmix(seed ^ base1) % w;
        while (base2 == base1);
        pw1.assign(1, 1), pw2.assign(1, 1);
    }
    static void ensure(std::size_t n) {  // pw1 の大きさを n + 1 以上に（基数は変えない）
        init();
        while (pw1.size() <= n) pw1.push_back(mul(pw1.back(), base1));
    }
    static void ensure2(std::size_t n) {  // pw2 も（2 次元用）
        ensure(n);
        while (pw2.size() <= n) pw2.push_back(mul(pw2.back(), base2));
    }
    static unsigned long long pow_mod(unsigned long long b, unsigned long long e) {
        unsigned long long r = 1;
        for (; e; e >>= 1, b = mul(b, b))
            if (e & 1) r = mul(r, b);
        return r;
    }
    static unsigned long long pow1(long long n) {  // base1^n（n が表の外なら繰り返し二乗）
        init();
        return (unsigned long long)n < pw1.size() ? pw1[(std::size_t)n] : pow_mod(base1, (unsigned long long)n);
    }
    static unsigned long long pow2(long long n) {
        init();
        return (unsigned long long)n < pw2.size() ? pw2[(std::size_t)n] : pow_mod(base2, (unsigned long long)n);
    }
    // 要素の値 → [1, MOD)。0〜2^60-1 はそのまま +1、それ以外は乱数で [2^60+1, MOD) に散らす
    static unsigned long long value(long long v) {
        if ((unsigned long long)v < (1ULL << 60)) return (unsigned long long)v + 1;
        init();
        return (1ULL << 60) + 1 + splitmix((unsigned long long)v ^ seed) % (MOD - (1ULL << 60) - 1);
    }
};
}  // namespace str_internal

// ---- 部分文字列のハッシュ値（長さも持つ。+ でつなげられる） ----
struct StrHash {
    unsigned long long value = 0;
    long long len = 0;
    friend StrHash operator+(const StrHash& a, const StrHash& b) {  // a のあとに b をつなげたもののハッシュ  O(1)（b.len が表の外なら O(log)）
        using C = str_internal::RHCore;
        return {C::add(C::mul(a.value, C::pow1(b.len)), b.value), a.len + b.len};
    }
    StrHash& operator+=(const StrHash& b) { return *this = *this + b; }
    StrHash repeat(long long k) const {  // k 回つなげたもの  O(log k)
        using C = str_internal::RHCore;
        if (k < 0 || (len > 0 && k > 4000000000000000000LL / len)) str_internal::fail("StrHash::repeat", "k ≥ 0 かつ長さ × k ≤ 4×10^18 で使う");
        // sum = 1 + P + … + P^(k-1)（P = base^len）を倍々で
        unsigned long long sum = 0, cur_sum = 1, cur_p = C::pow1(len);  // cur_sum = 1 + P + … + P^(2^b - 1)、cur_p = P^(2^b)
        for (long long e = k; e > 0; e >>= 1) {
            if (e & 1) sum = C::add(C::mul(sum, cur_p), cur_sum);
            cur_sum = C::add(cur_sum, C::mul(cur_sum, cur_p)), cur_p = C::mul(cur_p, cur_p);
        }
        return {C::mul(value, sum), len * k};
    }
    friend bool operator==(const StrHash& a, const StrHash& b) { return a.value == b.value && a.len == b.len; }
    friend bool operator!=(const StrHash& a, const StrHash& b) { return !(a == b); }
    friend bool operator<(const StrHash& a, const StrHash& b) { return a.len != b.len ? a.len < b.len : a.value < b.value; }  // map・sort 用（辞書順ではない）
    friend std::ostream& operator<<(std::ostream& os, const StrHash& h) { return os << "#" << h.value << "(len " << h.len << ")"; }  // デバッグ用
};
namespace std {
template <>
struct hash<StrHash> {  // unordered_map<StrHash, V> / unordered_set<StrHash> にそのまま使える
    size_t operator()(const StrHash& h) const noexcept { return (size_t)(h.value ^ ((unsigned long long)h.len * 0x9E3779B97F4A7C15ULL)); }
};
}  // namespace std

template <class Seq>
StrHash hash_of(const Seq& s0) {  // 列全体のハッシュ  O(N)
    using C = str_internal::RHCore;
    const auto& s = str_internal::seq(s0);
    C::init();
    unsigned long long h = 0;
    for (std::size_t i = 0; i < s.size(); ++i) h = C::add(C::mul(h, C::base1), C::value(str_internal::code_of(s[i])));
    return {h, (long long)s.size()};
}
template <class T>
StrHash hash_of_char(T c) { return {str_internal::RHCore::value(str_internal::code_of(c)), 1}; }  // 1 文字（1 要素）のハッシュ
inline StrHash empty_hash() { return {0, 0}; }

// ---- 列の部分文字列のハッシュを O(1) で ----
struct RollingHash {
    RollingHash() : h_(1, 0) { str_internal::RHCore::ensure(0); }  // 空の文字列と同じ
    template <class Seq>
    explicit RollingHash(const Seq& s0) {
        using C = str_internal::RHCore;
        const auto& s = str_internal::seq(s0);
        n_ = (int)s.size();
        C::ensure((std::size_t)n_);
        h_.resize(n_ + 1);
        h_[0] = 0;
        // 比べるときに使う元の値は数の列のときだけ持つ（char の値 0〜255 はハッシュから戻せる：1 文字のハッシュ = 値 + 1）
        constexpr bool narrow = std::is_same<str_internal::elem_t<Seq>, char>::value || std::is_same<str_internal::elem_t<Seq>, signed char>::value ||
                                std::is_same<str_internal::elem_t<Seq>, unsigned char>::value;
        if (!narrow) c64_.resize(n_);
        const unsigned long long b = C::base1;
        for (int i = 0; i < n_; ++i) {
            const long long c = str_internal::code_of(s[i]);
            if (!narrow) c64_[i] = c;
            h_[i + 1] = C::add(C::mul(h_[i], b), C::value(c));
        }
    }
    int size() const { return n_; }
    long long code(int i) const {  // i 文字目の値（char は 0〜255）
        if (i < 0 || i >= n_) str_internal::fail("RollingHash::code", "位置 " + std::to_string(i) + " は 0〜N-1 の外");
        return code_at(i);
    }
    StrHash get(int l, int r) const {  // s[l, r)  O(1)
        str_internal::check_range("RollingHash::get", l, r, n_);
        return {raw(l, r), (long long)(r - l)};
    }
    StrHash get_len(int l, int len) const { return get(l, l + len); }  // s[l, l + len)
    StrHash prefix(int r) const { return get(0, r); }                  // s[0, r)
    StrHash suffix(int l) const { return get(l, n_); }                 // s[l, N)
    StrHash all() const { return get(0, n_); }
    bool same(int l1, int r1, int l2, int r2) const { return get(l1, r1) == get(l2, r2); }  // s[l1, r1) == s[l2, r2) か  O(1)
    int lcp(int i, int j) const { return lcp(i, *this, j); }  // s[i..) と s[j..) の最長共通接頭辞の長さ  O(log N)
    int lcp(int i, const RollingHash& t, int j) const {       // s[i..) と t[j..) の最長共通接頭辞の長さ  O(log N)
        check_pos("RollingHash::lcp", i), t.check_pos("RollingHash::lcp", j);
        int lo = 0, hi = std::min(n_ - i, t.n_ - j);
        while (lo < hi) {
            const int mid = (lo + hi + 1) >> 1;
            if (raw(i, i + mid) == t.raw(j, j + mid)) lo = mid;
            else hi = mid - 1;
        }
        return lo;
    }
    int common_suffix(int i, int j) const { return common_suffix(i, *this, j); }  // s[0, i) と s[0, j) の最長共通接尾辞の長さ  O(log N)
    int common_suffix(int i, const RollingHash& t, int j) const {
        check_pos("RollingHash::common_suffix", i), t.check_pos("RollingHash::common_suffix", j);
        int lo = 0, hi = std::min(i, j);
        while (lo < hi) {
            const int mid = (lo + hi + 1) >> 1;
            if (raw(i - mid, i) == t.raw(j - mid, j)) lo = mid;
            else hi = mid - 1;
        }
        return lo;
    }
    int compare(int l1, int r1, int l2, int r2) const { return compare(l1, r1, *this, l2, r2); }  // s[l1, r1) と s[l2, r2) の辞書順（-1 / 0 / 1）  O(log N)
    int compare(int l1, int r1, const RollingHash& t, int l2, int r2) const {                      // s[l1, r1) と t[l2, r2)
        str_internal::check_range("RollingHash::compare", l1, r1, n_), str_internal::check_range("RollingHash::compare", l2, r2, t.n_);
        const int m = std::min(r1 - l1, r2 - l2);
        int lo = 0, hi = m;
        while (lo < hi) {
            const int mid = (lo + hi + 1) >> 1;
            if (raw(l1, l1 + mid) == t.raw(l2, l2 + mid)) lo = mid;
            else hi = mid - 1;
        }
        if (lo == m) return (r1 - l1) == (r2 - l2) ? 0 : ((r1 - l1) < (r2 - l2) ? -1 : 1);
        return code_at(l1 + lo) < t.code_at(l2 + lo) ? -1 : 1;
    }
    std::vector<int> find_all(const StrHash& p) const {  // ハッシュが p と同じ部分文字列の開始位置（昇順）  O(N)
        std::vector<int> r;
        if (p.len > n_) return r;
        const int m = (int)p.len;
        for (int i = 0; i + m <= n_; ++i)
            if (raw(i, i + m) == p.value) r.push_back(i);
        return r;
    }
    template <class Seq>
    std::vector<int> find_all(const Seq& pat) const { return find_all(hash_of(pat)); }  // pat が現れる開始位置（重なりも）
    template <class Seq>
    long long count(const Seq& pat) const { return (long long)find_all(pat).size(); }

   private:
    void check_pos(const char* func, int i) const {
        if (i < 0 || i > n_) str_internal::fail(func, "位置 " + std::to_string(i) + " は 0〜N（= " + std::to_string(n_) + "）の外");
    }
    unsigned long long raw(int l, int r) const {  // 範囲の確認なしのハッシュ値
        using C = str_internal::RHCore;
        return C::sub(h_[r], C::mul(h_[l], C::pw1[(std::size_t)(r - l)]));
    }
    long long code_at(int i) const { return c64_.empty() ? (long long)raw(i, i + 1) - 1 : c64_[i]; }
    int n_ = 0;
    std::vector<unsigned long long> h_;  // h_[i] = s[0, i) のハッシュ
    std::vector<long long> c64_;         // 元の値（数の列のとき。char の列なら空）
};

// ---- 2 次元のハッシュ（グリッドの長方形 [i1, i2) × [j1, j2) を O(1) で） ----
struct GridHash {
    unsigned long long value = 0;
    int h = 0, w = 0;
    friend bool operator==(const GridHash& a, const GridHash& b) { return a.value == b.value && a.h == b.h && a.w == b.w; }
    friend bool operator!=(const GridHash& a, const GridHash& b) { return !(a == b); }
    friend bool operator<(const GridHash& a, const GridHash& b) { return std::tie(a.h, a.w, a.value) < std::tie(b.h, b.w, b.value); }
};
namespace std {
template <>
struct hash<GridHash> {
    size_t operator()(const GridHash& g) const noexcept { return (size_t)(g.value ^ ((unsigned long long)g.h * 0x9E3779B97F4A7C15ULL) ^ ((unsigned long long)g.w << 32)); }
};
}  // namespace std
struct RollingHash2D {
    RollingHash2D() : P_(1, 0) { str_internal::RHCore::ensure2(0); }  // 0 × 0 のグリッドと同じ
    template <class Grid>
    explicit RollingHash2D(const Grid& g) {  // vector<string> でも vector<vector<int>> でも（行の長さはそろえる）
        using C = str_internal::RHCore;
        H_ = (int)g.size(), W_ = H_ ? (int)g[0].size() : 0;
        C::ensure2((std::size_t)std::max(H_, W_));
        P_.assign((std::size_t)(H_ + 1) * (W_ + 1), 0);
        const unsigned long long b12 = C::mul(C::base1, C::base2);
        for (int i = 0; i < H_; ++i) {
            if ((int)g[i].size() != W_) str_internal::fail("RollingHash2D", "行の長さがそろっていない");
            for (int j = 0; j < W_; ++j) {
                unsigned long long v = C::add(C::mul(at(i, j + 1), C::base2), C::mul(at(i + 1, j), C::base1));
                v = C::sub(v, C::mul(at(i, j), b12));
                at(i + 1, j + 1) = C::add(v, C::value(str_internal::code_of(g[i][j])));
            }
        }
    }
    int rows() const { return H_; }
    int cols() const { return W_; }
    GridHash get(int i1, int j1, int i2, int j2) const {  // 行 [i1, i2) × 列 [j1, j2)  O(1)
        str_internal::check_range("RollingHash2D::get", i1, i2, H_), str_internal::check_range("RollingHash2D::get", j1, j2, W_);
        using C = str_internal::RHCore;
        const unsigned long long ph = C::pw2[(std::size_t)(i2 - i1)], pw = C::pw1[(std::size_t)(j2 - j1)];
        unsigned long long v = C::sub(cat(i2, j2), C::mul(cat(i1, j2), ph));
        v = C::sub(v, C::mul(cat(i2, j1), pw));
        v = C::add(v, C::mul(C::mul(cat(i1, j1), ph), pw));
        return {v, i2 - i1, j2 - j1};
    }
    GridHash all() const { return get(0, 0, H_, W_); }
    template <class Grid>
    std::vector<std::pair<int, int>> find_all(const Grid& pat) const {  // pat と一致する左上の位置 (i, j) の一覧（行の順・列の順）  O(HW)
        std::vector<std::pair<int, int>> r;
        const RollingHash2D P(pat);
        const int h = P.rows(), w = P.cols();
        if (h > H_ || w > W_) return r;
        const GridHash want = P.all();
        for (int i = 0; i + h <= H_; ++i)
            for (int j = 0; j + w <= W_; ++j)
                if (get(i, j, i + h, j + w) == want) r.push_back({i, j});
        return r;
    }

   private:
    unsigned long long& at(int i, int j) { return P_[(std::size_t)i * (W_ + 1) + j]; }
    unsigned long long cat(int i, int j) const { return P_[(std::size_t)i * (W_ + 1) + j]; }
    int H_ = 0, W_ = 0;
    std::vector<unsigned long long> P_;  // P_[i][j] = 左上 [0, i) × [0, j) のハッシュ
};
#endif  // GG_ROLLING_HASH

#ifndef GG_STR_BASE
#define GG_STR_BASE
// ---- 文字列基盤（文字列第二世代）  文字の種類 Alphabet・文字の数と位置・ランレングス・変換・比較・区切り ----
// 列（string / vector<整数>）を受け取る関数はどちらでもよい。文字列リテラル（"abc"）も渡せる
// 範囲は半開区間 [l, r)（0 始まり）。位置が「ない」ときは -1。使える範囲の外の引数はメッセージを出して止める
namespace str_internal {
[[noreturn]] inline void fail(const char* func, const std::string& msg) {
    std::cerr << func << ": " << msg << "\n";
    std::abort();
}
// 文字列リテラルを std::string に（それ以外はそのまま参照を返す）
template <class T>
const T& seq(const T& s) { return s; }
template <std::size_t N>
std::string seq(const char (&s)[N]) { return std::string(s); }
inline std::string seq(const char* s) { return std::string(s); }
template <class T>
using seq_t = std::decay_t<decltype(seq(std::declval<const T&>()))>;
template <class T>
using elem_t = std::decay_t<decltype(std::declval<const seq_t<T>&>()[0])>;
// 要素を「値」に（char は unsigned char として 0〜255。負の添字にしない）
template <class T>
long long code_of(T x) {
    if constexpr (std::is_same<T, char>::value || std::is_same<T, signed char>::value || std::is_same<T, unsigned char>::value)
        return (long long)(unsigned char)x;
    else
        return (long long)x;
}
// 範囲の確認（止めるときの文字列作りは別の関数に分けて、確認そのものを軽くする）
[[noreturn]] __attribute__((noinline, cold)) inline void range_fail(const char* func, long long l, long long r, long long n) {
    fail(func, "範囲 [l, r) は 0 ≤ l ≤ r ≤ 長さ で指定する（l = " + std::to_string(l) + ", r = " + std::to_string(r) + ", 長さ = " + std::to_string(n) + "）");
}
inline void check_range(const char* func, long long l, long long r, long long n) {
    if (__builtin_expect(l < 0 || r > n || l > r, 0)) range_fail(func, l, r, n);
}
// a[lo, hi) の中で x 以上の最初の位置（a は昇順）。分岐の少ない二分探索
inline int lower_index(const int* a, int lo, int hi, int x) {
    int len = hi - lo;
    if (len <= 0) return lo;
    const int* b = a + lo;
    while (len > 1) {
        const int half = len >> 1;
        b += (b[half - 1] < x) ? half : 0;
        len -= half;
    }
    return (int)(b - a) + (*b < x);
}
// 昇順の a の中で x より大きい最初の位置（std::upper_bound と同じ。デバッグビルドでも O(log N)）
template <class T>
int upper_index(const std::vector<T>& a, const T& x) {
    int lo = 0, hi = (int)a.size();
    while (lo < hi) {
        const int mid = (int)((unsigned)(lo + hi) >> 1);
        if (!(x < a[mid])) lo = mid + 1;
        else hi = mid;
    }
    return lo;
}
// 数え上げの足し算：T が整数型なら溢れたら止める（mint・double などはそのまま）
template <class T>
T checked_add(const T& a, const T& b, const char* func) {
    if constexpr (std::is_integral<T>::value) {
        T r;
        if (__builtin_add_overflow(a, b, &r)) fail(func, "値が型に収まらない（答えを mod で求めるなら math_v2 の ModInt998 などを使う）");
        return r;
    } else {
        return a + b;
    }
}
// 数え上げの途中の足し算：T が整数型なら溢れたことを印 oa に残すだけ（途中で溢れても答えに使わなければ止めない。最後に印を見て止める）
template <class T>
void add_count(T& a, bool& oa, const T& b, bool ob) {
    if constexpr (std::is_integral<T>::value) {
        T r;
        oa = __builtin_add_overflow(a, b, &r) || oa || ob;
        a = r;
    } else {
        a = a + b;
    }
}
[[noreturn]] inline void count_overflow(const char* func) { fail(func, "値が型に収まらない（答えを mod で求めるなら math_v2 の ModInt998 などを使う）"); }
// 列 → 値の列（char は 0〜255）
template <class Seq>
std::vector<long long> codes_of(const Seq& s0) {
    const auto& s = seq(s0);
    std::vector<long long> c(s.size());
    for (std::size_t i = 0; i < s.size(); ++i) c[i] = code_of(s[i]);
    return c;
}
// 辞書順の比較（-1 / 0 / 1）。char は unsigned char として比べる（std::string と同じ順）
template <class A, class B>
int lex_compare(const A& a, const B& b) {
    const std::size_t n = a.size(), m = b.size(), k = n < m ? n : m;
    for (std::size_t i = 0; i < k; ++i) {
        const long long x = code_of(a[i]), y = code_of(b[i]);
        if (x != y) return x < y ? -1 : 1;
    }
    return n == m ? 0 : (n < m ? -1 : 1);
}
}  // namespace str_internal

// ---- 部分文字列の位置 [l, r) ----
struct StrSpan {
    int l = 0, r = 0;
    int len() const { return r - l; }
    bool empty() const { return r <= l; }
    friend bool operator==(const StrSpan& a, const StrSpan& b) { return a.l == b.l && a.r == b.r; }
    friend bool operator!=(const StrSpan& a, const StrSpan& b) { return !(a == b); }
    friend std::ostream& operator<<(std::ostream& os, const StrSpan& sp) { return os << "[" << sp.l << ", " << sp.r << ")"; }  // デバッグ用
};
template <class Seq>
str_internal::seq_t<Seq> substr_of(const Seq& s0, StrSpan sp) {  // s[sp.l, sp.r) のコピー（string なら string、vector なら vector）
    const auto& s = str_internal::seq(s0);
    str_internal::check_range("substr_of", sp.l, sp.r, (long long)s.size());
    return str_internal::seq_t<Seq>(s.begin() + sp.l, s.begin() + sp.r);
}

// ---- 文字の種類（Alphabet）：文字・値 → 番号 0〜size()-1。番号の順が辞書順として扱われる ----
//   Alphabet::lower()（a〜z。既定）/ upper()（A〜Z）/ digits()（0〜9）/ letters()（A〜Z, a〜z の 52 種。ASCII 順）
//   alnum()（0〜9, A〜Z, a〜z の 62 種。ASCII 順）/ ascii()（0〜127）/ bytes()（0〜255。日本語など 128 以上の char も）
//   of("()[]")（並べた文字だけ。並べた順が番号の順 = 辞書順）/ ints(K)（整数 0〜K-1）/ ints(lo, hi)（整数 [lo, hi)）
//   範囲外の文字・値を id に渡すとメッセージを出して止まる（contains で先に確かめられる）
struct Alphabet {
    Alphabet() : Alphabet(of("abcdefghijklmnopqrstuvwxyz")) {}
    static Alphabet lower() { return of("abcdefghijklmnopqrstuvwxyz"); }
    static Alphabet upper() { return of("ABCDEFGHIJKLMNOPQRSTUVWXYZ"); }
    static Alphabet digits() { return of("0123456789"); }
    static Alphabet letters() { return of("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"); }
    static Alphabet alnum() { return of("0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"); }
    static Alphabet ascii() { return code_range(0, 128); }
    static Alphabet bytes() { return code_range(0, 256); }
    static Alphabet of(const std::string& cs) {  // 並べた文字（同じ文字を 2 回並べると止まる）
        if (cs.empty()) str_internal::fail("Alphabet::of", "文字が 1 つもない");
        Alphabet A(0);
        A.kind_ = CHARS;
        for (std::size_t i = 0; i < cs.size(); ++i) {
            const unsigned char c = (unsigned char)cs[i];
            if (A.table_[c] >= 0) str_internal::fail("Alphabet::of", std::string("同じ文字が 2 回ある: '") + cs[i] + "'");
            A.table_[c] = (short)i;
            A.chars_[i] = c;
        }
        A.size_ = (int)cs.size();
        return A;
    }
    static Alphabet ints(long long K) { return ints(0, K); }  // 整数 0〜K-1
    static Alphabet ints(long long lo, long long hi) {         // 整数 [lo, hi)
        if (hi <= lo || (unsigned long long)hi - (unsigned long long)lo > 1000000000ULL)
            str_internal::fail("Alphabet::ints", "lo < hi かつ hi - lo ≤ 10^9 で指定する（大きな値は compress_values で 0〜K-1 に）");
        Alphabet A(0);
        A.kind_ = INTS;
        A.lo_ = lo;
        A.size_ = (int)(hi - lo);
        return A;
    }

    int size() const { return size_; }
    bool is_chars() const { return kind_ == CHARS; }  // 文字の表で決めたものか（of・lower など）。false なら ints
    template <class T>
    bool contains(T x) const { return index_of(str_internal::code_of(x)) >= 0; }
    template <class T>
    int id(T x) const {  // 文字・値の番号（範囲外は止める）
        const long long v = str_internal::code_of(x);
        const int r = index_of(v);
        if (r < 0) out_of_range(v);
        return r;
    }
    long long value(int i) const {  // 番号 i の値（文字なら文字コード 0〜255）
        check_index(i);
        return kind_ == CHARS ? (long long)chars_[i] : lo_ + i;
    }
    char ch(int i) const {  // 番号 i の文字（ints で 0〜255 の外の値なら止める）
        const long long v = value(i);
        if (v < 0 || v > 255) str_internal::fail("Alphabet::ch", "番号 " + std::to_string(i) + " の値 " + std::to_string(v) + " は文字にできない");
        return (char)(unsigned char)v;
    }
    std::string chars() const {  // 文字を番号の順に並べたもの（ints なら 0〜255 の範囲のときだけ）
        std::string s;
        for (int i = 0; i < size_; ++i) s += ch(i);
        return s;
    }

   private:
    enum Kind : unsigned char { CHARS, INTS };
    explicit Alphabet(int) { table_.fill(-1); }
    static Alphabet code_range(int lo, int hi) {
        Alphabet A(0);
        A.kind_ = CHARS;
        for (int c = lo; c < hi; ++c) A.table_[c] = (short)(c - lo), A.chars_[c - lo] = (unsigned char)c;
        A.size_ = hi - lo;
        return A;
    }
    int index_of(long long v) const {
        if (kind_ == CHARS) return (unsigned long long)v < 256ULL ? table_[(std::size_t)v] : -1;
        const unsigned long long d = (unsigned long long)v - (unsigned long long)lo_;
        return d < (unsigned long long)size_ ? (int)d : -1;
    }
    void check_index(int i) const {
        if (i < 0 || i >= size_) str_internal::fail("Alphabet", "番号 " + std::to_string(i) + " は 0〜" + std::to_string(size_ - 1) + " の外");
    }
    [[noreturn]] void out_of_range(long long v) const {
        std::string what = "値 " + std::to_string(v) + " は";
        if (kind_ == CHARS && v >= 32 && v < 127) what = std::string("文字 '") + (char)v + "'（コード " + std::to_string(v) + "）は";
        str_internal::fail("Alphabet::id", what + "文字の種類に入っていない（種類の数 " + std::to_string(size_) + "）");
    }
    Kind kind_ = CHARS;
    int size_ = 0;
    long long lo_ = 0;
    std::array<short, 256> table_{};
    std::array<unsigned char, 256> chars_{};
};

template <class Seq>
std::vector<int> to_ids(const Seq& s0, const Alphabet& A = Alphabet::lower()) {  // 各文字の番号の列
    const auto& s = str_internal::seq(s0);
    std::vector<int> r(s.size());
    for (std::size_t i = 0; i < s.size(); ++i) r[i] = A.id(s[i]);
    return r;
}
inline std::string from_ids(const std::vector<int>& ids, const Alphabet& A = Alphabet::lower()) {  // 番号の列 → 文字列
    std::string s(ids.size(), '\0');
    for (std::size_t i = 0; i < ids.size(); ++i) s[i] = A.ch(ids[i]);
    return s;
}

// ---- 数の列を 0〜K-1 に圧縮（大小の順を保つ）：Alphabet::ints(K) と組み合わせる ----
struct CompressedValues {
    std::vector<int> ids;            // 各要素の番号（0〜K-1。値の小さい順）
    std::vector<long long> values;   // 番号 → 元の値（昇順。K = values.size()）
    int size() const { return (int)values.size(); }
    Alphabet alphabet() const { return Alphabet::ints(values.empty() ? 1 : (long long)values.size()); }
};
template <class Seq>
CompressedValues compress_values(const Seq& s0) {
    const auto& s = str_internal::seq(s0);
    const std::size_t n = s.size();
    std::vector<std::pair<long long, int>> vi(n);
    for (std::size_t i = 0; i < n; ++i) vi[i] = {str_internal::code_of(s[i]), (int)i};
    std::sort(vi.begin(), vi.end());
    CompressedValues c;
    c.ids.resize(n);
    for (std::size_t i = 0; i < n; ++i) {
        if (i == 0 || vi[i].first != vi[i - 1].first) c.values.push_back(vi[i].first);
        c.ids[vi[i].second] = (int)c.values.size() - 1;
    }
    return c;
}

// ---- 文字の数 ----
template <class Seq>
std::vector<int> count_chars(const Seq& s0, const Alphabet& A = Alphabet::lower()) {  // 番号ごとの個数（長さ A.size()）
    const auto& s = str_internal::seq(s0);
    std::vector<int> c(A.size(), 0);
    for (std::size_t i = 0; i < s.size(); ++i) ++c[A.id(s[i])];
    return c;
}

// ---- 区間の文字の数・次 / 前の出現を O(1) で（前処理 O(Nσ)・メモリ (N+1)σ 個の int。N = 2×10^5, σ = 26 で約 21MB） ----
//   メモリを抑えたいとき・σ が大きいときは CharPositions（O(log N)・メモリ O(N + σ)）
struct CharPrefix {
    CharPrefix() = default;
    template <class Seq>
    explicit CharPrefix(const Seq& s0, const Alphabet& A = Alphabet::lower()) : A_(A) {
        const auto& s = str_internal::seq(s0);
        n_ = (int)s.size(), sg_ = A_.size();
        if ((long long)(n_ + 1) * sg_ > 600000000LL) str_internal::fail("CharPrefix", "(N + 1)σ が 6×10^8 を超える（メモリが足りない。CharPositions を使う）");
        cnt_.assign((std::size_t)(n_ + 1) * sg_, 0);
        std::vector<int> id(n_);
        start_.assign(sg_ + 1, 0);
        for (int i = 0; i < n_; ++i) {
            const int* src = cnt_.data() + (std::size_t)i * sg_;
            int* dst = cnt_.data() + (std::size_t)(i + 1) * sg_;
            std::copy(src, src + sg_, dst);
            id[i] = A_.id(s[i]);
            ++dst[id[i]], ++start_[id[i] + 1];
        }
        for (int k = 0; k < sg_; ++k) start_[k + 1] += start_[k];
        pos_.resize(n_);
        std::vector<int> cur(start_.begin(), start_.end() - 1);
        for (int i = 0; i < n_; ++i) pos_[cur[id[i]]++] = i;
    }
    int size() const { return n_; }
    const Alphabet& alphabet() const { return A_; }
    template <class C>
    int count(int l, int r, C c) const {  // [l, r) にある c の数  O(1)
        str_internal::check_range("CharPrefix::count", l, r, n_);
        const int k = A_.id(c);
        return cnt_[(std::size_t)r * sg_ + k] - cnt_[(std::size_t)l * sg_ + k];
    }
    template <class C>
    bool has(int l, int r, C c) const { return count(l, r, c) > 0; }  // [l, r) に c があるか  O(1)
    std::vector<int> counts(int l, int r) const {  // [l, r) の番号ごとの数  O(σ)
        str_internal::check_range("CharPrefix::counts", l, r, n_);
        std::vector<int> c(sg_);
        for (int k = 0; k < sg_; ++k) c[k] = cnt_[(std::size_t)r * sg_ + k] - cnt_[(std::size_t)l * sg_ + k];
        return c;
    }
    int distinct(int l, int r) const {  // [l, r) にある文字の種類の数  O(σ)
        str_internal::check_range("CharPrefix::distinct", l, r, n_);
        int d = 0;
        for (int k = 0; k < sg_; ++k) d += cnt_[(std::size_t)r * sg_ + k] != cnt_[(std::size_t)l * sg_ + k];
        return d;
    }
    template <class C>
    int next(int i, C c) const {  // i 以上で最初の c の位置（なければ -1。0 ≤ i ≤ N）  O(1)
        check_pos("CharPrefix::next", i);
        const int k = A_.id(c), m = cnt_[(std::size_t)i * sg_ + k];
        return start_[k] + m < start_[k + 1] ? pos_[start_[k] + m] : -1;
    }
    template <class C>
    int prev(int i, C c) const {  // i より前（[0, i)）で最後の c の位置（なければ -1。0 ≤ i ≤ N）  O(1)
        check_pos("CharPrefix::prev", i);
        const int k = A_.id(c), m = cnt_[(std::size_t)i * sg_ + k];
        return m > 0 ? pos_[start_[k] + m - 1] : -1;
    }
    template <class C>
    int kth(C c, long long k) const {  // c の k 番目（0 から）の位置（なければ -1）  O(1)
        const int id = A_.id(c);
        return (k >= 0 && k < start_[id + 1] - start_[id]) ? pos_[start_[id] + k] : -1;
    }
    template <class C>
    int total(C c) const {  // 全体での c の数  O(1)
        const int k = A_.id(c);
        return start_[k + 1] - start_[k];
    }

   private:
    void check_pos(const char* func, int i) const {
        if (i < 0 || i > n_) str_internal::fail(func, "位置 " + std::to_string(i) + " は 0〜N（= " + std::to_string(n_) + "）の外");
    }
    Alphabet A_;
    int n_ = 0, sg_ = 0;
    std::vector<int> cnt_;           // cnt_[i * σ + k] = [0, i) にある番号 k の数
    std::vector<int> start_, pos_;   // 番号ごとの出現位置（pos_[start_[k] .. start_[k+1])）
};

// ---- 文字の出現位置（前処理 O(N + σ)・メモリ O(N + σ)。σ が大きくても使える） ----
struct CharPositions {
    CharPositions() = default;
    template <class Seq>
    explicit CharPositions(const Seq& s0, const Alphabet& A = Alphabet::lower()) : A_(A) {
        const auto& s = str_internal::seq(s0);
        n_ = (int)s.size();
        const int sg = A_.size();
        std::vector<int> id(n_);
        start_.assign(sg + 1, 0);
        for (int i = 0; i < n_; ++i) id[i] = A_.id(s[i]), ++start_[id[i] + 1];
        for (int k = 0; k < sg; ++k) start_[k + 1] += start_[k];
        pos_.resize(n_), occ_.resize(n_);
        std::vector<int> cur(start_.begin(), start_.end() - 1);
        for (int i = 0; i < n_; ++i) occ_[i] = cur[id[i]] - start_[id[i]], pos_[cur[id[i]]++] = i;
    }
    int size() const { return n_; }
    const Alphabet& alphabet() const { return A_; }
    template <class C>
    int next(int i, C c) const {  // i 以上で最初の c の位置（なければ -1。0 ≤ i ≤ N）  O(log N)
        check_pos("CharPositions::next", i);
        const int k = A_.id(c), j = str_internal::lower_index(pos_.data(), start_[k], start_[k + 1], i);
        return j < start_[k + 1] ? pos_[j] : -1;
    }
    template <class C>
    int prev(int i, C c) const {  // i より前（[0, i)）で最後の c の位置（なければ -1。0 ≤ i ≤ N）  O(log N)
        check_pos("CharPositions::prev", i);
        const int k = A_.id(c), j = str_internal::lower_index(pos_.data(), start_[k], start_[k + 1], i);
        return j > start_[k] ? pos_[j - 1] : -1;
    }
    template <class C>
    int kth(C c, long long k) const {  // c の k 番目（0 から）の位置（なければ -1）  O(1)
        const int id = A_.id(c);
        return (k >= 0 && k < start_[id + 1] - start_[id]) ? pos_[start_[id] + k] : -1;
    }
    template <class C>
    int count(int l, int r, C c) const {  // [l, r) にある c の数  O(log N)
        str_internal::check_range("CharPositions::count", l, r, n_);
        const int k = A_.id(c);
        return str_internal::lower_index(pos_.data(), start_[k], start_[k + 1], r) - str_internal::lower_index(pos_.data(), start_[k], start_[k + 1], l);
    }
    template <class C>
    bool has(int l, int r, C c) const { return count(l, r, c) > 0; }
    template <class C>
    int total(C c) const {  // 全体での c の数  O(1)
        const int k = A_.id(c);
        return start_[k + 1] - start_[k];
    }
    template <class C>
    std::vector<int> positions(C c) const {  // c の位置の一覧（昇順）
        const int k = A_.id(c);
        return std::vector<int>(pos_.begin() + start_[k], pos_.begin() + start_[k + 1]);
    }
    int occurrence_index(int i) const {  // s[i] が「その文字の何番目（0 から）の出現」か  O(1)
        if (i < 0 || i >= n_) str_internal::fail("CharPositions::occurrence_index", "位置 " + std::to_string(i) + " は 0〜N-1 の外");
        return occ_[i];
    }

   private:
    void check_pos(const char* func, int i) const {
        if (i < 0 || i > n_) str_internal::fail(func, "位置 " + std::to_string(i) + " は 0〜N（= " + std::to_string(n_) + "）の外");
    }
    Alphabet A_;
    int n_ = 0;
    std::vector<int> start_, pos_, occ_;
};

// ---- ランレングス圧縮：連続する同じ値を {value, len, start} にまとめる ----
template <class T>
struct RunBlock {
    T value{};
    long long len = 0;    // 続く長さ
    long long start = 0;  // 元の列での開始位置
    friend std::ostream& operator<<(std::ostream& os, const RunBlock& b) { return os << "{" << b.value << ", " << b.len << ", " << b.start << "}"; }  // デバッグ用
};
template <class Seq>
std::vector<RunBlock<str_internal::elem_t<Seq>>> run_length(const Seq& s0) {
    const auto& s = str_internal::seq(s0);
    std::vector<RunBlock<str_internal::elem_t<Seq>>> r;
    std::size_t cnt = s.empty() ? 0 : 1;
    for (std::size_t i = 1; i < s.size(); ++i) cnt += !(s[i] == s[i - 1]);
    r.reserve(cnt);
    for (std::size_t i = 0, n = s.size(); i < n;) {
        std::size_t j = i + 1;
        while (j < n && s[j] == s[i]) ++j;
        r.push_back({s[i], (long long)(j - i), (long long)i});
        i = j;
    }
    return r;
}
template <class T>
auto run_length_decode(const std::vector<RunBlock<T>>& runs) {  // 元に戻す（T = char なら string、それ以外は vector<T>）
    using R = std::conditional_t<std::is_same<T, char>::value, std::string, std::vector<T>>;
    long long total = 0;
    for (const auto& b : runs) {
        if (b.len < 0) str_internal::fail("run_length_decode", "len が負");
        total += b.len;
        if (total > 2000000000LL) str_internal::fail("run_length_decode", "元の長さが 2×10^9 を超える");
    }
    R r;
    r.reserve((std::size_t)total);
    for (const auto& b : runs) r.insert(r.end(), (std::size_t)b.len, b.value);
    return r;
}

// ---- 変換 ----
inline std::string to_lower(std::string s) {  // 英大文字を小文字に（ほかはそのまま）
    for (char& c : s)
        if ('A' <= c && c <= 'Z') c = (char)(c - 'A' + 'a');
    return s;
}
inline std::string to_upper(std::string s) {
    for (char& c : s)
        if ('a' <= c && c <= 'z') c = (char)(c - 'a' + 'A');
    return s;
}
inline std::string swap_case(std::string s) {  // 大文字と小文字を入れ替える
    for (char& c : s) {
        if ('A' <= c && c <= 'Z') c = (char)(c - 'A' + 'a');
        else if ('a' <= c && c <= 'z') c = (char)(c - 'a' + 'A');
    }
    return s;
}
template <class Seq>
str_internal::seq_t<Seq> caesar_shift(const Seq& s0, long long k, const Alphabet& A = Alphabet::lower()) {  // A の文字を番号で k ずらす（負も。A にない文字はそのまま）
    auto s = str_internal::seq_t<Seq>(str_internal::seq(s0));
    const long long sg = A.size(), sh = ((k % sg) + sg) % sg;
    for (auto& c : s)
        if (A.contains(c)) c = (std::decay_t<decltype(c)>)A.value((int)((A.id(c) + sh) % sg));
    return s;
}
// 文字の置き換えを重ねる：replace(c, d) で「今 c になっている文字を全部 d に」。apply(s) で s に当てはめる（ABC342 C）
struct CharMap {
    explicit CharMap(const Alphabet& A = Alphabet::lower()) : A_(A), f_(A.size()) {
        for (int i = 0; i < A_.size(); ++i) f_[i] = i;
    }
    template <class C, class D>
    void replace(C c, D d) {  // 今 c になっているものを全部 d に  O(σ)
        const int x = A_.id(c), y = A_.id(d);
        if (x == y) return;
        for (int& v : f_)
            if (v == x) v = y;
    }
    template <class C, class D>
    void set(C c, D d) { f_[A_.id(c)] = A_.id(d); }  // 元の文字 c の行き先を d にする（ほかは変えない）
    template <class C>
    auto operator()(C c) const {  // 元の文字 c の行き先（A にない文字はそのまま）
        return A_.contains(c) ? (C)A_.value(f_[A_.id(c)]) : c;
    }
    template <class Seq>
    str_internal::seq_t<Seq> apply(const Seq& s0) const {  // s の各文字を行き先に（A にない文字はそのまま）
        auto s = str_internal::seq_t<Seq>(str_internal::seq(s0));
        for (auto& c : s) c = (*this)(c);
        return s;
    }
    CharMap compose(const CharMap& g) const {  // 「この置き換えのあとに g」を 1 つにしたもの
        CharMap h(A_);
        for (int i = 0; i < A_.size(); ++i) h.f_[i] = g.f_[f_[i]];
        return h;
    }
    const std::vector<int>& table() const { return f_; }  // 番号 → 行き先の番号

   private:
    Alphabet A_;
    std::vector<int> f_;
};
template <class Seq>
str_internal::seq_t<Seq> reversed(const Seq& s0) {  // 反転したコピー
    const auto& s = str_internal::seq(s0);
    return str_internal::seq_t<Seq>(s.rbegin(), s.rend());
}
template <class Seq>
str_internal::seq_t<Seq> min_with_reverse(const Seq& s0) {  // s と反転のうち辞書順で小さい方（向きを区別しないときの代表。ABC310 C）
    auto s = str_internal::seq_t<Seq>(str_internal::seq(s0));
    auto t = reversed(s);
    return str_internal::lex_compare(t, s) < 0 ? t : s;
}
template <class Seq>
str_internal::seq_t<Seq> sorted_chars(const Seq& s0) {  // 並べ替えたもの（並べ替えて同じかの比較用）
    auto s = str_internal::seq_t<Seq>(str_internal::seq(s0));
    std::sort(s.begin(), s.end(), [](const auto& a, const auto& b) { return str_internal::code_of(a) < str_internal::code_of(b); });
    return s;
}
template <class Seq>
str_internal::seq_t<Seq> repeat_string(const Seq& s0, long long k) {  // s を k 回つなげる
    const auto& s = str_internal::seq(s0);
    if (k < 0) str_internal::fail("repeat_string", "k ≥ 0 で使う");
    if (!s.empty() && k > 2000000000LL / (long long)s.size()) str_internal::fail("repeat_string", "長さが 2×10^9 を超える");
    str_internal::seq_t<Seq> r;
    r.reserve((std::size_t)(k * (long long)s.size()));
    for (long long i = 0; i < k; ++i) r.insert(r.end(), s.begin(), s.end());
    return r;
}

// ---- 比較 ----
template <class A, class B>
int common_prefix_length(const A& a0, const B& b0) {  // 共通の接頭辞の長さ
    const auto& a = str_internal::seq(a0);
    const auto& b = str_internal::seq(b0);
    const std::size_t k = a.size() < b.size() ? a.size() : b.size();
    std::size_t i = 0;
    while (i < k && str_internal::code_of(a[i]) == str_internal::code_of(b[i])) ++i;
    return (int)i;
}
template <class A, class B>
int common_suffix_length(const A& a0, const B& b0) {  // 共通の接尾辞の長さ
    const auto& a = str_internal::seq(a0);
    const auto& b = str_internal::seq(b0);
    const std::size_t k = a.size() < b.size() ? a.size() : b.size();
    std::size_t i = 0;
    while (i < k && str_internal::code_of(a[a.size() - 1 - i]) == str_internal::code_of(b[b.size() - 1 - i])) ++i;
    return (int)i;
}
template <class A, class B>
int hamming_distance(const A& a0, const B& b0) {  // 同じ位置で違う文字の数（長さが同じこと）
    const auto& a = str_internal::seq(a0);
    const auto& b = str_internal::seq(b0);
    if (a.size() != b.size()) str_internal::fail("hamming_distance", "長さが違う");
    int d = 0;
    for (std::size_t i = 0; i < a.size(); ++i) d += str_internal::code_of(a[i]) != str_internal::code_of(b[i]);
    return d;
}
template <class T, class S>
bool is_subsequence(const T& t0, const S& s0) {  // t が s の部分列（飛び飛びでよい）か  O(|s|)。何度も聞くなら 06_subsequence の NextTable
    const auto& t = str_internal::seq(t0);
    const auto& s = str_internal::seq(s0);
    std::size_t j = 0;
    for (std::size_t i = 0; i < s.size() && j < t.size(); ++i)
        if (str_internal::code_of(s[i]) == str_internal::code_of(t[j])) ++j;
    return j == t.size();
}
template <class Seq>
bool is_palindrome(const Seq& s0) {  // 全体が回文か  O(N)。区間を何度も聞くなら 04_palindrome の Palindromes
    const auto& s = str_internal::seq(s0);
    for (std::size_t i = 0, n = s.size(); i < n / 2; ++i)
        if (!(s[i] == s[n - 1 - i])) return false;
    return true;
}

// ---- 区切り・つなげる ----
inline std::vector<std::string> split(const std::string& s, char sep, bool skip_empty = false) {  // sep で区切る（空の部分も残す。skip_empty で除く）
    std::vector<std::string> r;
    std::size_t b = 0;
    for (std::size_t i = 0; i <= s.size(); ++i)
        if (i == s.size() || s[i] == sep) {
            if (!(skip_empty && i == b)) r.push_back(s.substr(b, i - b));
            b = i + 1;
        }
    return r;
}
template <class V>
std::string join(const std::vector<V>& vs, const std::string& sep = " ") {  // 区切りを挟んでつなげる（数も文字列にする。答えの出力に）
    std::string r;
    for (std::size_t i = 0; i < vs.size(); ++i) {
        if (i) r += sep;
        if constexpr (std::is_same<V, std::string>::value) r += vs[i];
        else if constexpr (std::is_same<V, char>::value) r += vs[i];
        else if constexpr (std::is_integral<V>::value) r += std::to_string(vs[i]);
        else {
            std::ostringstream os;
            os << vs[i];
            r += os.str();
        }
    }
    return r;
}
#endif  // GG_STR_BASE
#ifndef GG_PATTERN_MATCH
#define GG_PATTERN_MATCH
// ---- パターン検索（文字列第二世代）  Z・KMP・出現・周期・境界・回転・重なり・KMP オートマトン ----
// text・pat は string でも vector<整数> でもよい（混ぜてもよい。値で比べる）。文字列リテラルも渡せる
namespace str_internal {
// 列（string・vector<整数>）をそのまま読む（コピーしない）。要素は code_of で比べる
template <class S>
std::vector<int> z_codes(const S& s) {
    const int n = (int)s.size();
    std::vector<int> z(n);
    if (n == 0) return z;
    z[0] = n;
    for (int i = 1, l = 0, r = 0; i < n; ++i) {
        int k = i < r ? std::min(r - i, z[i - l]) : 0;
        while (i + k < n && s[k] == s[i + k]) ++k;
        z[i] = k;
        if (i + k > r) l = i, r = i + k;
    }
    return z;
}
template <class S>
std::vector<int> pi_codes(const S& s) {
    const int n = (int)s.size();
    std::vector<int> pi(n);
    for (int i = 1; i < n; ++i) {
        int k = pi[i - 1];
        while (k > 0 && !(s[i] == s[k])) k = pi[k - 1];
        pi[i] = k + (s[i] == s[k]);
    }
    return pi;
}
// text[from..) の中の pat の出現開始位置（text での位置。重なりも）。pat が空なら from〜|text|
template <class T, class P>
std::vector<int> kmp_find(const T& t, const P& p, int from = 0) {
    const int n = (int)t.size(), m = (int)p.size();
    std::vector<int> r;
    if (m == 0) {
        for (int i = from; i <= n; ++i) r.push_back(i);
        return r;
    }
    const std::vector<int> pi = pi_codes(p);
    for (int i = from, k = 0; i < n; ++i) {
        const long long c = code_of(t[i]);
        while (k > 0 && (k == m || c != code_of(p[k]))) k = pi[k - 1];
        if (c == code_of(p[k])) ++k;
        if (k == m) r.push_back(i - m + 1);
    }
    return r;
}
}  // namespace str_internal

template <class Seq>
std::vector<int> z_array(const Seq& s) { return str_internal::z_codes(str_internal::seq(s)); }  // z[i] = s と s[i..) の最長共通接頭辞（z[0] = |s|）  O(N)
template <class Seq>
std::vector<int> prefix_function(const Seq& s) { return str_internal::pi_codes(str_internal::seq(s)); }  // π[i] = s[0, i] の最長の真の境界の長さ  O(N)
template <class T, class P>
std::vector<int> match_lengths(const T& text, const P& pat) {  // 各 i で text[i..) と pat の最長共通接頭辞の長さ（長さ |text|）  O(N + M)
    const auto& t = str_internal::seq(text);
    const auto& p = str_internal::seq(pat);
    const int n = (int)t.size(), m = (int)p.size();
    const std::vector<int> zp = str_internal::z_codes(p);
    std::vector<int> res(n);
    for (int i = 0, l = 0, r = 0; i < n; ++i) {  // [l, r) は text[l, r) == pat[0, r - l) となる右端が最大の区間
        int k = 0;
        if (i < r) k = std::min(r - i, zp[i - l]);
        if (i >= r || zp[i - l] >= r - i) {
            while (k < m && i + k < n && str_internal::code_of(t[i + k]) == str_internal::code_of(p[k])) ++k;
        }
        res[i] = k;
        if (i + k > r) l = i, r = i + k;
    }
    return res;
}
template <class T, class P>
std::vector<int> find_all(const T& text, const P& pat) {  // pat が現れる開始位置（重なりも。昇順。pat が空なら 0〜|text|）  O(N + M)
    return str_internal::kmp_find(str_internal::seq(text), str_internal::seq(pat));
}
template <class T, class P>
long long count_occurrences(const T& text, const P& pat) { return (long long)find_all(text, pat).size(); }  // 重なりも数える
template <class T, class P>
int find_first(const T& text, const P& pat, int from = 0) {  // from 以上で最初に現れる位置（なければ -1）  O(N + M)
    const auto& t = str_internal::seq(text);
    if (from < 0 || from > (int)t.size()) str_internal::fail("find_first", "from は 0〜|text| で指定する");
    const auto r = str_internal::kmp_find(t, str_internal::seq(pat), from);
    return r.empty() ? -1 : r[0];
}
template <class T, class P>
int find_last(const T& text, const P& pat) {  // 最後に現れる位置（なければ -1）
    const auto r = find_all(text, pat);
    return r.empty() ? -1 : r.back();
}
template <class T, class P>
long long count_non_overlapping(const T& text, const P& pat) {  // 左から重ならないように取った数（pat は空でない）
    const auto& p = str_internal::seq(pat);
    if (p.empty()) str_internal::fail("count_non_overlapping", "pat が空");
    long long c = 0, last = -1;
    for (int i : str_internal::kmp_find(str_internal::seq(text), p))
        if (i >= last) ++c, last = i + (long long)p.size();
    return c;
}
template <class T, class F, class U>
str_internal::seq_t<T> replace_all(const T& text0, const F& from, const U& to0) {  // from を左から重ならないように全部 to に（from は空でない）  O(N + M)
    const auto& text = str_internal::seq(text0);
    const auto& to = str_internal::seq(to0);
    const auto& f = str_internal::seq(from);
    if (f.empty()) str_internal::fail("replace_all", "from が空");
    str_internal::seq_t<T> r;
    std::size_t pos = 0;
    for (int i : str_internal::kmp_find(text, f)) {
        if ((std::size_t)i < pos) continue;
        r.insert(r.end(), text.begin() + pos, text.begin() + i);
        for (std::size_t j = 0; j < to.size(); ++j) r.push_back((typename str_internal::seq_t<T>::value_type)str_internal::code_of(to[j]));
        pos = (std::size_t)i + f.size();
    }
    r.insert(r.end(), text.begin() + pos, text.end());
    return r;
}

// ---- 周期・境界 ----
template <class Seq>
std::vector<int> all_borders(const Seq& s) {  // 0 < k < |s| で s[0, k) == s[|s|-k, |s|) となる k（大きい順）  O(N)
    const auto pi = prefix_function(s);
    std::vector<int> r;
    if (pi.empty()) return r;
    for (int k = pi.back(); k > 0; k = pi[k - 1]) r.push_back(k);
    return r;
}
template <class Seq>
int min_period(const Seq& s) {  // 最小の周期 p（s[i] == s[i + p]。割り切れなくてよい。空なら 0）  O(N)
    const auto pi = prefix_function(s);
    return pi.empty() ? 0 : (int)pi.size() - pi.back();
}
template <class Seq>
std::vector<int> all_periods(const Seq& s) {  // 1〜|s| のすべての周期（小さい順。|s| 自身も含む）  O(N)
    const int n = (int)str_internal::seq(s).size();
    std::vector<int> r;
    for (int k : all_borders(s)) r.push_back(n - k);
    if (n > 0) r.push_back(n);
    return r;
}
template <class Seq>
int repetition_root(const Seq& s) {  // s = t^k（t を k 回）となる最短の t の長さ（空なら 0）  O(N)
    const int n = (int)str_internal::seq(s).size(), p = min_period(s);
    return (n > 0 && n % p == 0) ? p : n;
}
template <class Seq>
std::vector<int> prefix_min_periods(const Seq& s) {  // r[i] = s[0, i) の最小の周期（長さ |s| + 1。r[0] = 0）  O(N)
    const auto pi = prefix_function(s);
    std::vector<int> r(pi.size() + 1, 0);
    for (std::size_t i = 0; i < pi.size(); ++i) r[i + 1] = (int)(i + 1) - pi[i];
    return r;
}

// ---- 回転・重なり ----
template <class A, class B>
std::vector<int> rotation_offsets(const A& a, const B& b) {  // a を左に k 回まわす（a[k..] + a[..k]）と b になる k（0 ≤ k < |a|。昇順。a, b がともに空なら {0}）  O(N)
    const auto x = str_internal::codes_of(a), y = str_internal::codes_of(b);
    std::vector<int> r;
    if (x.size() != y.size()) return r;
    if (x.empty()) return {0};
    std::vector<long long> xx(x);
    xx.insert(xx.end(), x.begin(), x.end() - 1);
    return str_internal::kmp_find(xx, y);
}
template <class A, class B>
bool is_rotation(const A& a, const B& b) { return !rotation_offsets(a, b).empty(); }  // 回転して一致するか
template <class T, class P>
std::vector<int> find_all_cyclic(const T& text, const P& pat) {  // text を無限に繰り返した列で pat が始まる位置 i（0 ≤ i < |text|。|pat| > |text| でもよい）  O(N + M)
    const auto t = str_internal::codes_of(text), p = str_internal::codes_of(pat);
    const std::size_t n = t.size(), m = p.size();
    if (n == 0) return {};
    std::vector<long long> tt;
    tt.reserve(n + m);
    while (tt.size() < n + (m > 0 ? m - 1 : 0)) tt.insert(tt.end(), t.begin(), t.end());
    std::vector<int> r;
    for (int i : str_internal::kmp_find(tt, p))
        if ((std::size_t)i < n) r.push_back(i);
    return r;
}
template <class A, class B>
int max_overlap(const A& a, const B& b) {  // a の接尾辞と b の接頭辞が一致する最大の長さ（a + b を短くつなげるとき）  O(N + M)
    const auto x = str_internal::codes_of(a), y = str_internal::codes_of(b);
    const int m = (int)y.size();
    if (m == 0) return 0;
    const std::vector<int> pi = str_internal::pi_codes(y);
    int k = 0;
    for (long long c : x) {
        if (k == m) k = pi[k - 1];
        while (k > 0 && c != y[k]) k = pi[k - 1];
        if (c == y[k]) ++k;
    }
    return k;
}

// ---- KMP オートマトン：状態 = pat の先頭何文字まで一致しているか（0〜|pat|）。|pat| に着いたら pat が現れた ----
struct KmpAutomaton {
    KmpAutomaton() = default;
    template <class P>
    explicit KmpAutomaton(const P& pat0, const Alphabet& A = Alphabet::lower()) : A_(A) {
        const auto& pat = str_internal::seq(pat0);
        m_ = (int)pat.size();
        sg_ = A_.size();
        std::vector<int> p(m_);
        for (int i = 0; i < m_; ++i) p[i] = A_.id(pat[i]);
        nx_.assign((std::size_t)(m_ + 1) * sg_, 0);
        int fail = 0;  // 状態 q の「1 つ短い一致」
        for (int q = 0; q <= m_; ++q) {
            for (int c = 0; c < sg_; ++c) nx_[(std::size_t)q * sg_ + c] = (q < m_ && p[q] == c) ? q + 1 : (q == 0 ? 0 : nx_[(std::size_t)fail * sg_ + c]);
            if (q > 0 && q < m_) fail = nx_[(std::size_t)fail * sg_ + p[q]];
        }
    }
    int size() const { return m_ + 1; }  // 状態の数（0〜|pat|）
    template <class C>
    int next(int state, C c) const { return next_id(state, A_.id(c)); }  // 状態 state で文字 c を読んだあとの状態  O(1)
    int next_id(int state, int id) const {
        if (state < 0 || state > m_ || id < 0 || id >= sg_) str_internal::fail("KmpAutomaton::next", "状態か文字の番号が範囲外");
        return nx_[(std::size_t)state * sg_ + id];
    }
    bool accepting(int state) const { return state == m_; }  // pat がちょうど現れたところか
    const Alphabet& alphabet() const { return A_; }

   private:
    Alphabet A_;
    int m_ = 0, sg_ = 0;
    std::vector<int> nx_;
};
#endif  // GG_PATTERN_MATCH

#ifndef GG_STR_BASE
#define GG_STR_BASE
// ---- 文字列基盤（文字列第二世代）  文字の種類 Alphabet・文字の数と位置・ランレングス・変換・比較・区切り ----
// 列（string / vector<整数>）を受け取る関数はどちらでもよい。文字列リテラル（"abc"）も渡せる
// 範囲は半開区間 [l, r)（0 始まり）。位置が「ない」ときは -1。使える範囲の外の引数はメッセージを出して止める
namespace str_internal {
[[noreturn]] inline void fail(const char* func, const std::string& msg) {
    std::cerr << func << ": " << msg << "\n";
    std::abort();
}
// 文字列リテラルを std::string に（それ以外はそのまま参照を返す）
template <class T>
const T& seq(const T& s) { return s; }
template <std::size_t N>
std::string seq(const char (&s)[N]) { return std::string(s); }
inline std::string seq(const char* s) { return std::string(s); }
template <class T>
using seq_t = std::decay_t<decltype(seq(std::declval<const T&>()))>;
template <class T>
using elem_t = std::decay_t<decltype(std::declval<const seq_t<T>&>()[0])>;
// 要素を「値」に（char は unsigned char として 0〜255。負の添字にしない）
template <class T>
long long code_of(T x) {
    if constexpr (std::is_same<T, char>::value || std::is_same<T, signed char>::value || std::is_same<T, unsigned char>::value)
        return (long long)(unsigned char)x;
    else
        return (long long)x;
}
// 範囲の確認（止めるときの文字列作りは別の関数に分けて、確認そのものを軽くする）
[[noreturn]] __attribute__((noinline, cold)) inline void range_fail(const char* func, long long l, long long r, long long n) {
    fail(func, "範囲 [l, r) は 0 ≤ l ≤ r ≤ 長さ で指定する（l = " + std::to_string(l) + ", r = " + std::to_string(r) + ", 長さ = " + std::to_string(n) + "）");
}
inline void check_range(const char* func, long long l, long long r, long long n) {
    if (__builtin_expect(l < 0 || r > n || l > r, 0)) range_fail(func, l, r, n);
}
// a[lo, hi) の中で x 以上の最初の位置（a は昇順）。分岐の少ない二分探索
inline int lower_index(const int* a, int lo, int hi, int x) {
    int len = hi - lo;
    if (len <= 0) return lo;
    const int* b = a + lo;
    while (len > 1) {
        const int half = len >> 1;
        b += (b[half - 1] < x) ? half : 0;
        len -= half;
    }
    return (int)(b - a) + (*b < x);
}
// 昇順の a の中で x より大きい最初の位置（std::upper_bound と同じ。デバッグビルドでも O(log N)）
template <class T>
int upper_index(const std::vector<T>& a, const T& x) {
    int lo = 0, hi = (int)a.size();
    while (lo < hi) {
        const int mid = (int)((unsigned)(lo + hi) >> 1);
        if (!(x < a[mid])) lo = mid + 1;
        else hi = mid;
    }
    return lo;
}
// 数え上げの足し算：T が整数型なら溢れたら止める（mint・double などはそのまま）
template <class T>
T checked_add(const T& a, const T& b, const char* func) {
    if constexpr (std::is_integral<T>::value) {
        T r;
        if (__builtin_add_overflow(a, b, &r)) fail(func, "値が型に収まらない（答えを mod で求めるなら math_v2 の ModInt998 などを使う）");
        return r;
    } else {
        return a + b;
    }
}
// 数え上げの途中の足し算：T が整数型なら溢れたことを印 oa に残すだけ（途中で溢れても答えに使わなければ止めない。最後に印を見て止める）
template <class T>
void add_count(T& a, bool& oa, const T& b, bool ob) {
    if constexpr (std::is_integral<T>::value) {
        T r;
        oa = __builtin_add_overflow(a, b, &r) || oa || ob;
        a = r;
    } else {
        a = a + b;
    }
}
[[noreturn]] inline void count_overflow(const char* func) { fail(func, "値が型に収まらない（答えを mod で求めるなら math_v2 の ModInt998 などを使う）"); }
// 列 → 値の列（char は 0〜255）
template <class Seq>
std::vector<long long> codes_of(const Seq& s0) {
    const auto& s = seq(s0);
    std::vector<long long> c(s.size());
    for (std::size_t i = 0; i < s.size(); ++i) c[i] = code_of(s[i]);
    return c;
}
// 辞書順の比較（-1 / 0 / 1）。char は unsigned char として比べる（std::string と同じ順）
template <class A, class B>
int lex_compare(const A& a, const B& b) {
    const std::size_t n = a.size(), m = b.size(), k = n < m ? n : m;
    for (std::size_t i = 0; i < k; ++i) {
        const long long x = code_of(a[i]), y = code_of(b[i]);
        if (x != y) return x < y ? -1 : 1;
    }
    return n == m ? 0 : (n < m ? -1 : 1);
}
}  // namespace str_internal

// ---- 部分文字列の位置 [l, r) ----
struct StrSpan {
    int l = 0, r = 0;
    int len() const { return r - l; }
    bool empty() const { return r <= l; }
    friend bool operator==(const StrSpan& a, const StrSpan& b) { return a.l == b.l && a.r == b.r; }
    friend bool operator!=(const StrSpan& a, const StrSpan& b) { return !(a == b); }
    friend std::ostream& operator<<(std::ostream& os, const StrSpan& sp) { return os << "[" << sp.l << ", " << sp.r << ")"; }  // デバッグ用
};
template <class Seq>
str_internal::seq_t<Seq> substr_of(const Seq& s0, StrSpan sp) {  // s[sp.l, sp.r) のコピー（string なら string、vector なら vector）
    const auto& s = str_internal::seq(s0);
    str_internal::check_range("substr_of", sp.l, sp.r, (long long)s.size());
    return str_internal::seq_t<Seq>(s.begin() + sp.l, s.begin() + sp.r);
}

// ---- 文字の種類（Alphabet）：文字・値 → 番号 0〜size()-1。番号の順が辞書順として扱われる ----
//   Alphabet::lower()（a〜z。既定）/ upper()（A〜Z）/ digits()（0〜9）/ letters()（A〜Z, a〜z の 52 種。ASCII 順）
//   alnum()（0〜9, A〜Z, a〜z の 62 種。ASCII 順）/ ascii()（0〜127）/ bytes()（0〜255。日本語など 128 以上の char も）
//   of("()[]")（並べた文字だけ。並べた順が番号の順 = 辞書順）/ ints(K)（整数 0〜K-1）/ ints(lo, hi)（整数 [lo, hi)）
//   範囲外の文字・値を id に渡すとメッセージを出して止まる（contains で先に確かめられる）
struct Alphabet {
    Alphabet() : Alphabet(of("abcdefghijklmnopqrstuvwxyz")) {}
    static Alphabet lower() { return of("abcdefghijklmnopqrstuvwxyz"); }
    static Alphabet upper() { return of("ABCDEFGHIJKLMNOPQRSTUVWXYZ"); }
    static Alphabet digits() { return of("0123456789"); }
    static Alphabet letters() { return of("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"); }
    static Alphabet alnum() { return of("0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"); }
    static Alphabet ascii() { return code_range(0, 128); }
    static Alphabet bytes() { return code_range(0, 256); }
    static Alphabet of(const std::string& cs) {  // 並べた文字（同じ文字を 2 回並べると止まる）
        if (cs.empty()) str_internal::fail("Alphabet::of", "文字が 1 つもない");
        Alphabet A(0);
        A.kind_ = CHARS;
        for (std::size_t i = 0; i < cs.size(); ++i) {
            const unsigned char c = (unsigned char)cs[i];
            if (A.table_[c] >= 0) str_internal::fail("Alphabet::of", std::string("同じ文字が 2 回ある: '") + cs[i] + "'");
            A.table_[c] = (short)i;
            A.chars_[i] = c;
        }
        A.size_ = (int)cs.size();
        return A;
    }
    static Alphabet ints(long long K) { return ints(0, K); }  // 整数 0〜K-1
    static Alphabet ints(long long lo, long long hi) {         // 整数 [lo, hi)
        if (hi <= lo || (unsigned long long)hi - (unsigned long long)lo > 1000000000ULL)
            str_internal::fail("Alphabet::ints", "lo < hi かつ hi - lo ≤ 10^9 で指定する（大きな値は compress_values で 0〜K-1 に）");
        Alphabet A(0);
        A.kind_ = INTS;
        A.lo_ = lo;
        A.size_ = (int)(hi - lo);
        return A;
    }

    int size() const { return size_; }
    bool is_chars() const { return kind_ == CHARS; }  // 文字の表で決めたものか（of・lower など）。false なら ints
    template <class T>
    bool contains(T x) const { return index_of(str_internal::code_of(x)) >= 0; }
    template <class T>
    int id(T x) const {  // 文字・値の番号（範囲外は止める）
        const long long v = str_internal::code_of(x);
        const int r = index_of(v);
        if (r < 0) out_of_range(v);
        return r;
    }
    long long value(int i) const {  // 番号 i の値（文字なら文字コード 0〜255）
        check_index(i);
        return kind_ == CHARS ? (long long)chars_[i] : lo_ + i;
    }
    char ch(int i) const {  // 番号 i の文字（ints で 0〜255 の外の値なら止める）
        const long long v = value(i);
        if (v < 0 || v > 255) str_internal::fail("Alphabet::ch", "番号 " + std::to_string(i) + " の値 " + std::to_string(v) + " は文字にできない");
        return (char)(unsigned char)v;
    }
    std::string chars() const {  // 文字を番号の順に並べたもの（ints なら 0〜255 の範囲のときだけ）
        std::string s;
        for (int i = 0; i < size_; ++i) s += ch(i);
        return s;
    }

   private:
    enum Kind : unsigned char { CHARS, INTS };
    explicit Alphabet(int) { table_.fill(-1); }
    static Alphabet code_range(int lo, int hi) {
        Alphabet A(0);
        A.kind_ = CHARS;
        for (int c = lo; c < hi; ++c) A.table_[c] = (short)(c - lo), A.chars_[c - lo] = (unsigned char)c;
        A.size_ = hi - lo;
        return A;
    }
    int index_of(long long v) const {
        if (kind_ == CHARS) return (unsigned long long)v < 256ULL ? table_[(std::size_t)v] : -1;
        const unsigned long long d = (unsigned long long)v - (unsigned long long)lo_;
        return d < (unsigned long long)size_ ? (int)d : -1;
    }
    void check_index(int i) const {
        if (i < 0 || i >= size_) str_internal::fail("Alphabet", "番号 " + std::to_string(i) + " は 0〜" + std::to_string(size_ - 1) + " の外");
    }
    [[noreturn]] void out_of_range(long long v) const {
        std::string what = "値 " + std::to_string(v) + " は";
        if (kind_ == CHARS && v >= 32 && v < 127) what = std::string("文字 '") + (char)v + "'（コード " + std::to_string(v) + "）は";
        str_internal::fail("Alphabet::id", what + "文字の種類に入っていない（種類の数 " + std::to_string(size_) + "）");
    }
    Kind kind_ = CHARS;
    int size_ = 0;
    long long lo_ = 0;
    std::array<short, 256> table_{};
    std::array<unsigned char, 256> chars_{};
};

template <class Seq>
std::vector<int> to_ids(const Seq& s0, const Alphabet& A = Alphabet::lower()) {  // 各文字の番号の列
    const auto& s = str_internal::seq(s0);
    std::vector<int> r(s.size());
    for (std::size_t i = 0; i < s.size(); ++i) r[i] = A.id(s[i]);
    return r;
}
inline std::string from_ids(const std::vector<int>& ids, const Alphabet& A = Alphabet::lower()) {  // 番号の列 → 文字列
    std::string s(ids.size(), '\0');
    for (std::size_t i = 0; i < ids.size(); ++i) s[i] = A.ch(ids[i]);
    return s;
}

// ---- 数の列を 0〜K-1 に圧縮（大小の順を保つ）：Alphabet::ints(K) と組み合わせる ----
struct CompressedValues {
    std::vector<int> ids;            // 各要素の番号（0〜K-1。値の小さい順）
    std::vector<long long> values;   // 番号 → 元の値（昇順。K = values.size()）
    int size() const { return (int)values.size(); }
    Alphabet alphabet() const { return Alphabet::ints(values.empty() ? 1 : (long long)values.size()); }
};
template <class Seq>
CompressedValues compress_values(const Seq& s0) {
    const auto& s = str_internal::seq(s0);
    const std::size_t n = s.size();
    std::vector<std::pair<long long, int>> vi(n);
    for (std::size_t i = 0; i < n; ++i) vi[i] = {str_internal::code_of(s[i]), (int)i};
    std::sort(vi.begin(), vi.end());
    CompressedValues c;
    c.ids.resize(n);
    for (std::size_t i = 0; i < n; ++i) {
        if (i == 0 || vi[i].first != vi[i - 1].first) c.values.push_back(vi[i].first);
        c.ids[vi[i].second] = (int)c.values.size() - 1;
    }
    return c;
}

// ---- 文字の数 ----
template <class Seq>
std::vector<int> count_chars(const Seq& s0, const Alphabet& A = Alphabet::lower()) {  // 番号ごとの個数（長さ A.size()）
    const auto& s = str_internal::seq(s0);
    std::vector<int> c(A.size(), 0);
    for (std::size_t i = 0; i < s.size(); ++i) ++c[A.id(s[i])];
    return c;
}

// ---- 区間の文字の数・次 / 前の出現を O(1) で（前処理 O(Nσ)・メモリ (N+1)σ 個の int。N = 2×10^5, σ = 26 で約 21MB） ----
//   メモリを抑えたいとき・σ が大きいときは CharPositions（O(log N)・メモリ O(N + σ)）
struct CharPrefix {
    CharPrefix() = default;
    template <class Seq>
    explicit CharPrefix(const Seq& s0, const Alphabet& A = Alphabet::lower()) : A_(A) {
        const auto& s = str_internal::seq(s0);
        n_ = (int)s.size(), sg_ = A_.size();
        if ((long long)(n_ + 1) * sg_ > 600000000LL) str_internal::fail("CharPrefix", "(N + 1)σ が 6×10^8 を超える（メモリが足りない。CharPositions を使う）");
        cnt_.assign((std::size_t)(n_ + 1) * sg_, 0);
        std::vector<int> id(n_);
        start_.assign(sg_ + 1, 0);
        for (int i = 0; i < n_; ++i) {
            const int* src = cnt_.data() + (std::size_t)i * sg_;
            int* dst = cnt_.data() + (std::size_t)(i + 1) * sg_;
            std::copy(src, src + sg_, dst);
            id[i] = A_.id(s[i]);
            ++dst[id[i]], ++start_[id[i] + 1];
        }
        for (int k = 0; k < sg_; ++k) start_[k + 1] += start_[k];
        pos_.resize(n_);
        std::vector<int> cur(start_.begin(), start_.end() - 1);
        for (int i = 0; i < n_; ++i) pos_[cur[id[i]]++] = i;
    }
    int size() const { return n_; }
    const Alphabet& alphabet() const { return A_; }
    template <class C>
    int count(int l, int r, C c) const {  // [l, r) にある c の数  O(1)
        str_internal::check_range("CharPrefix::count", l, r, n_);
        const int k = A_.id(c);
        return cnt_[(std::size_t)r * sg_ + k] - cnt_[(std::size_t)l * sg_ + k];
    }
    template <class C>
    bool has(int l, int r, C c) const { return count(l, r, c) > 0; }  // [l, r) に c があるか  O(1)
    std::vector<int> counts(int l, int r) const {  // [l, r) の番号ごとの数  O(σ)
        str_internal::check_range("CharPrefix::counts", l, r, n_);
        std::vector<int> c(sg_);
        for (int k = 0; k < sg_; ++k) c[k] = cnt_[(std::size_t)r * sg_ + k] - cnt_[(std::size_t)l * sg_ + k];
        return c;
    }
    int distinct(int l, int r) const {  // [l, r) にある文字の種類の数  O(σ)
        str_internal::check_range("CharPrefix::distinct", l, r, n_);
        int d = 0;
        for (int k = 0; k < sg_; ++k) d += cnt_[(std::size_t)r * sg_ + k] != cnt_[(std::size_t)l * sg_ + k];
        return d;
    }
    template <class C>
    int next(int i, C c) const {  // i 以上で最初の c の位置（なければ -1。0 ≤ i ≤ N）  O(1)
        check_pos("CharPrefix::next", i);
        const int k = A_.id(c), m = cnt_[(std::size_t)i * sg_ + k];
        return start_[k] + m < start_[k + 1] ? pos_[start_[k] + m] : -1;
    }
    template <class C>
    int prev(int i, C c) const {  // i より前（[0, i)）で最後の c の位置（なければ -1。0 ≤ i ≤ N）  O(1)
        check_pos("CharPrefix::prev", i);
        const int k = A_.id(c), m = cnt_[(std::size_t)i * sg_ + k];
        return m > 0 ? pos_[start_[k] + m - 1] : -1;
    }
    template <class C>
    int kth(C c, long long k) const {  // c の k 番目（0 から）の位置（なければ -1）  O(1)
        const int id = A_.id(c);
        return (k >= 0 && k < start_[id + 1] - start_[id]) ? pos_[start_[id] + k] : -1;
    }
    template <class C>
    int total(C c) const {  // 全体での c の数  O(1)
        const int k = A_.id(c);
        return start_[k + 1] - start_[k];
    }

   private:
    void check_pos(const char* func, int i) const {
        if (i < 0 || i > n_) str_internal::fail(func, "位置 " + std::to_string(i) + " は 0〜N（= " + std::to_string(n_) + "）の外");
    }
    Alphabet A_;
    int n_ = 0, sg_ = 0;
    std::vector<int> cnt_;           // cnt_[i * σ + k] = [0, i) にある番号 k の数
    std::vector<int> start_, pos_;   // 番号ごとの出現位置（pos_[start_[k] .. start_[k+1])）
};

// ---- 文字の出現位置（前処理 O(N + σ)・メモリ O(N + σ)。σ が大きくても使える） ----
struct CharPositions {
    CharPositions() = default;
    template <class Seq>
    explicit CharPositions(const Seq& s0, const Alphabet& A = Alphabet::lower()) : A_(A) {
        const auto& s = str_internal::seq(s0);
        n_ = (int)s.size();
        const int sg = A_.size();
        std::vector<int> id(n_);
        start_.assign(sg + 1, 0);
        for (int i = 0; i < n_; ++i) id[i] = A_.id(s[i]), ++start_[id[i] + 1];
        for (int k = 0; k < sg; ++k) start_[k + 1] += start_[k];
        pos_.resize(n_), occ_.resize(n_);
        std::vector<int> cur(start_.begin(), start_.end() - 1);
        for (int i = 0; i < n_; ++i) occ_[i] = cur[id[i]] - start_[id[i]], pos_[cur[id[i]]++] = i;
    }
    int size() const { return n_; }
    const Alphabet& alphabet() const { return A_; }
    template <class C>
    int next(int i, C c) const {  // i 以上で最初の c の位置（なければ -1。0 ≤ i ≤ N）  O(log N)
        check_pos("CharPositions::next", i);
        const int k = A_.id(c), j = str_internal::lower_index(pos_.data(), start_[k], start_[k + 1], i);
        return j < start_[k + 1] ? pos_[j] : -1;
    }
    template <class C>
    int prev(int i, C c) const {  // i より前（[0, i)）で最後の c の位置（なければ -1。0 ≤ i ≤ N）  O(log N)
        check_pos("CharPositions::prev", i);
        const int k = A_.id(c), j = str_internal::lower_index(pos_.data(), start_[k], start_[k + 1], i);
        return j > start_[k] ? pos_[j - 1] : -1;
    }
    template <class C>
    int kth(C c, long long k) const {  // c の k 番目（0 から）の位置（なければ -1）  O(1)
        const int id = A_.id(c);
        return (k >= 0 && k < start_[id + 1] - start_[id]) ? pos_[start_[id] + k] : -1;
    }
    template <class C>
    int count(int l, int r, C c) const {  // [l, r) にある c の数  O(log N)
        str_internal::check_range("CharPositions::count", l, r, n_);
        const int k = A_.id(c);
        return str_internal::lower_index(pos_.data(), start_[k], start_[k + 1], r) - str_internal::lower_index(pos_.data(), start_[k], start_[k + 1], l);
    }
    template <class C>
    bool has(int l, int r, C c) const { return count(l, r, c) > 0; }
    template <class C>
    int total(C c) const {  // 全体での c の数  O(1)
        const int k = A_.id(c);
        return start_[k + 1] - start_[k];
    }
    template <class C>
    std::vector<int> positions(C c) const {  // c の位置の一覧（昇順）
        const int k = A_.id(c);
        return std::vector<int>(pos_.begin() + start_[k], pos_.begin() + start_[k + 1]);
    }
    int occurrence_index(int i) const {  // s[i] が「その文字の何番目（0 から）の出現」か  O(1)
        if (i < 0 || i >= n_) str_internal::fail("CharPositions::occurrence_index", "位置 " + std::to_string(i) + " は 0〜N-1 の外");
        return occ_[i];
    }

   private:
    void check_pos(const char* func, int i) const {
        if (i < 0 || i > n_) str_internal::fail(func, "位置 " + std::to_string(i) + " は 0〜N（= " + std::to_string(n_) + "）の外");
    }
    Alphabet A_;
    int n_ = 0;
    std::vector<int> start_, pos_, occ_;
};

// ---- ランレングス圧縮：連続する同じ値を {value, len, start} にまとめる ----
template <class T>
struct RunBlock {
    T value{};
    long long len = 0;    // 続く長さ
    long long start = 0;  // 元の列での開始位置
    friend std::ostream& operator<<(std::ostream& os, const RunBlock& b) { return os << "{" << b.value << ", " << b.len << ", " << b.start << "}"; }  // デバッグ用
};
template <class Seq>
std::vector<RunBlock<str_internal::elem_t<Seq>>> run_length(const Seq& s0) {
    const auto& s = str_internal::seq(s0);
    std::vector<RunBlock<str_internal::elem_t<Seq>>> r;
    std::size_t cnt = s.empty() ? 0 : 1;
    for (std::size_t i = 1; i < s.size(); ++i) cnt += !(s[i] == s[i - 1]);
    r.reserve(cnt);
    for (std::size_t i = 0, n = s.size(); i < n;) {
        std::size_t j = i + 1;
        while (j < n && s[j] == s[i]) ++j;
        r.push_back({s[i], (long long)(j - i), (long long)i});
        i = j;
    }
    return r;
}
template <class T>
auto run_length_decode(const std::vector<RunBlock<T>>& runs) {  // 元に戻す（T = char なら string、それ以外は vector<T>）
    using R = std::conditional_t<std::is_same<T, char>::value, std::string, std::vector<T>>;
    long long total = 0;
    for (const auto& b : runs) {
        if (b.len < 0) str_internal::fail("run_length_decode", "len が負");
        total += b.len;
        if (total > 2000000000LL) str_internal::fail("run_length_decode", "元の長さが 2×10^9 を超える");
    }
    R r;
    r.reserve((std::size_t)total);
    for (const auto& b : runs) r.insert(r.end(), (std::size_t)b.len, b.value);
    return r;
}

// ---- 変換 ----
inline std::string to_lower(std::string s) {  // 英大文字を小文字に（ほかはそのまま）
    for (char& c : s)
        if ('A' <= c && c <= 'Z') c = (char)(c - 'A' + 'a');
    return s;
}
inline std::string to_upper(std::string s) {
    for (char& c : s)
        if ('a' <= c && c <= 'z') c = (char)(c - 'a' + 'A');
    return s;
}
inline std::string swap_case(std::string s) {  // 大文字と小文字を入れ替える
    for (char& c : s) {
        if ('A' <= c && c <= 'Z') c = (char)(c - 'A' + 'a');
        else if ('a' <= c && c <= 'z') c = (char)(c - 'a' + 'A');
    }
    return s;
}
template <class Seq>
str_internal::seq_t<Seq> caesar_shift(const Seq& s0, long long k, const Alphabet& A = Alphabet::lower()) {  // A の文字を番号で k ずらす（負も。A にない文字はそのまま）
    auto s = str_internal::seq_t<Seq>(str_internal::seq(s0));
    const long long sg = A.size(), sh = ((k % sg) + sg) % sg;
    for (auto& c : s)
        if (A.contains(c)) c = (std::decay_t<decltype(c)>)A.value((int)((A.id(c) + sh) % sg));
    return s;
}
// 文字の置き換えを重ねる：replace(c, d) で「今 c になっている文字を全部 d に」。apply(s) で s に当てはめる（ABC342 C）
struct CharMap {
    explicit CharMap(const Alphabet& A = Alphabet::lower()) : A_(A), f_(A.size()) {
        for (int i = 0; i < A_.size(); ++i) f_[i] = i;
    }
    template <class C, class D>
    void replace(C c, D d) {  // 今 c になっているものを全部 d に  O(σ)
        const int x = A_.id(c), y = A_.id(d);
        if (x == y) return;
        for (int& v : f_)
            if (v == x) v = y;
    }
    template <class C, class D>
    void set(C c, D d) { f_[A_.id(c)] = A_.id(d); }  // 元の文字 c の行き先を d にする（ほかは変えない）
    template <class C>
    auto operator()(C c) const {  // 元の文字 c の行き先（A にない文字はそのまま）
        return A_.contains(c) ? (C)A_.value(f_[A_.id(c)]) : c;
    }
    template <class Seq>
    str_internal::seq_t<Seq> apply(const Seq& s0) const {  // s の各文字を行き先に（A にない文字はそのまま）
        auto s = str_internal::seq_t<Seq>(str_internal::seq(s0));
        for (auto& c : s) c = (*this)(c);
        return s;
    }
    CharMap compose(const CharMap& g) const {  // 「この置き換えのあとに g」を 1 つにしたもの
        CharMap h(A_);
        for (int i = 0; i < A_.size(); ++i) h.f_[i] = g.f_[f_[i]];
        return h;
    }
    const std::vector<int>& table() const { return f_; }  // 番号 → 行き先の番号

   private:
    Alphabet A_;
    std::vector<int> f_;
};
template <class Seq>
str_internal::seq_t<Seq> reversed(const Seq& s0) {  // 反転したコピー
    const auto& s = str_internal::seq(s0);
    return str_internal::seq_t<Seq>(s.rbegin(), s.rend());
}
template <class Seq>
str_internal::seq_t<Seq> min_with_reverse(const Seq& s0) {  // s と反転のうち辞書順で小さい方（向きを区別しないときの代表。ABC310 C）
    auto s = str_internal::seq_t<Seq>(str_internal::seq(s0));
    auto t = reversed(s);
    return str_internal::lex_compare(t, s) < 0 ? t : s;
}
template <class Seq>
str_internal::seq_t<Seq> sorted_chars(const Seq& s0) {  // 並べ替えたもの（並べ替えて同じかの比較用）
    auto s = str_internal::seq_t<Seq>(str_internal::seq(s0));
    std::sort(s.begin(), s.end(), [](const auto& a, const auto& b) { return str_internal::code_of(a) < str_internal::code_of(b); });
    return s;
}
template <class Seq>
str_internal::seq_t<Seq> repeat_string(const Seq& s0, long long k) {  // s を k 回つなげる
    const auto& s = str_internal::seq(s0);
    if (k < 0) str_internal::fail("repeat_string", "k ≥ 0 で使う");
    if (!s.empty() && k > 2000000000LL / (long long)s.size()) str_internal::fail("repeat_string", "長さが 2×10^9 を超える");
    str_internal::seq_t<Seq> r;
    r.reserve((std::size_t)(k * (long long)s.size()));
    for (long long i = 0; i < k; ++i) r.insert(r.end(), s.begin(), s.end());
    return r;
}

// ---- 比較 ----
template <class A, class B>
int common_prefix_length(const A& a0, const B& b0) {  // 共通の接頭辞の長さ
    const auto& a = str_internal::seq(a0);
    const auto& b = str_internal::seq(b0);
    const std::size_t k = a.size() < b.size() ? a.size() : b.size();
    std::size_t i = 0;
    while (i < k && str_internal::code_of(a[i]) == str_internal::code_of(b[i])) ++i;
    return (int)i;
}
template <class A, class B>
int common_suffix_length(const A& a0, const B& b0) {  // 共通の接尾辞の長さ
    const auto& a = str_internal::seq(a0);
    const auto& b = str_internal::seq(b0);
    const std::size_t k = a.size() < b.size() ? a.size() : b.size();
    std::size_t i = 0;
    while (i < k && str_internal::code_of(a[a.size() - 1 - i]) == str_internal::code_of(b[b.size() - 1 - i])) ++i;
    return (int)i;
}
template <class A, class B>
int hamming_distance(const A& a0, const B& b0) {  // 同じ位置で違う文字の数（長さが同じこと）
    const auto& a = str_internal::seq(a0);
    const auto& b = str_internal::seq(b0);
    if (a.size() != b.size()) str_internal::fail("hamming_distance", "長さが違う");
    int d = 0;
    for (std::size_t i = 0; i < a.size(); ++i) d += str_internal::code_of(a[i]) != str_internal::code_of(b[i]);
    return d;
}
template <class T, class S>
bool is_subsequence(const T& t0, const S& s0) {  // t が s の部分列（飛び飛びでよい）か  O(|s|)。何度も聞くなら 06_subsequence の NextTable
    const auto& t = str_internal::seq(t0);
    const auto& s = str_internal::seq(s0);
    std::size_t j = 0;
    for (std::size_t i = 0; i < s.size() && j < t.size(); ++i)
        if (str_internal::code_of(s[i]) == str_internal::code_of(t[j])) ++j;
    return j == t.size();
}
template <class Seq>
bool is_palindrome(const Seq& s0) {  // 全体が回文か  O(N)。区間を何度も聞くなら 04_palindrome の Palindromes
    const auto& s = str_internal::seq(s0);
    for (std::size_t i = 0, n = s.size(); i < n / 2; ++i)
        if (!(s[i] == s[n - 1 - i])) return false;
    return true;
}

// ---- 区切り・つなげる ----
inline std::vector<std::string> split(const std::string& s, char sep, bool skip_empty = false) {  // sep で区切る（空の部分も残す。skip_empty で除く）
    std::vector<std::string> r;
    std::size_t b = 0;
    for (std::size_t i = 0; i <= s.size(); ++i)
        if (i == s.size() || s[i] == sep) {
            if (!(skip_empty && i == b)) r.push_back(s.substr(b, i - b));
            b = i + 1;
        }
    return r;
}
template <class V>
std::string join(const std::vector<V>& vs, const std::string& sep = " ") {  // 区切りを挟んでつなげる（数も文字列にする。答えの出力に）
    std::string r;
    for (std::size_t i = 0; i < vs.size(); ++i) {
        if (i) r += sep;
        if constexpr (std::is_same<V, std::string>::value) r += vs[i];
        else if constexpr (std::is_same<V, char>::value) r += vs[i];
        else if constexpr (std::is_integral<V>::value) r += std::to_string(vs[i]);
        else {
            std::ostringstream os;
            os << vs[i];
            r += os.str();
        }
    }
    return r;
}
#endif  // GG_STR_BASE
#ifndef GG_SUFFIX_ARRAY
#define GG_SUFFIX_ARRAY
// ---- 接尾辞配列（文字列第二世代）  SA-IS・LCP（Kasai）・O(1) の LCP と比較・出現・異なる部分文字列・k 番目・複数の文字列 ----
// ACL の suffix_array / lcp_array と名前を変えている（make_suffix_array / make_lcp_array。using namespace atcoder; と一緒でも衝突しない）
namespace str_internal {
inline std::vector<int> sa_naive(const std::vector<int>& s) {
    const int n = (int)s.size();
    std::vector<int> sa(n);
    std::iota(sa.begin(), sa.end(), 0);
    std::sort(sa.begin(), sa.end(), [&](int l, int r) {
        if (l == r) return false;
        while (l < n && r < n) {
            if (s[l] != s[r]) return s[l] < s[r];
            ++l, ++r;
        }
        return l == n;
    });
    return sa;
}
inline std::vector<int> sa_doubling(const std::vector<int>& s) {
    const int n = (int)s.size();
    std::vector<int> sa(n), rnk = s, tmp(n);
    std::iota(sa.begin(), sa.end(), 0);
    for (int k = 1; k < n; k *= 2) {
        auto cmp = [&](int x, int y) {
            if (rnk[x] != rnk[y]) return rnk[x] < rnk[y];
            const int rx = x + k < n ? rnk[x + k] : -1, ry = y + k < n ? rnk[y + k] : -1;
            return rx < ry;
        };
        std::sort(sa.begin(), sa.end(), cmp);
        tmp[sa[0]] = 0;
        for (int i = 1; i < n; ++i) tmp[sa[i]] = tmp[sa[i - 1]] + (cmp(sa[i - 1], sa[i]) ? 1 : 0);
        std::swap(tmp, rnk);
    }
    return sa;
}
// SA-IS（s の値は 0〜upper）
inline std::vector<int> sa_is(const std::vector<int>& s, int upper) {
    const int n = (int)s.size();
    if (n == 0) return {};
    if (n == 1) return {0};
    if (n == 2) return s[0] < s[1] ? std::vector<int>{0, 1} : std::vector<int>{1, 0};
    if (n < 10) return sa_naive(s);
    if (n < 40) return sa_doubling(s);
    std::vector<int> sa(n);
    std::vector<bool> ls(n);  // true = S 型（vector<bool> の方が速い：induce で飛び飛びに読むので小さい方がよい）
    for (int i = n - 2; i >= 0; --i) ls[i] = (s[i] == s[i + 1]) ? ls[i + 1] : (s[i] < s[i + 1]);
    std::vector<int> sum_l(upper + 1), sum_s(upper + 1);
    for (int i = 0; i < n; ++i) {
        if (!ls[i]) ++sum_s[s[i]];
        else ++sum_l[s[i] + 1];
    }
    for (int i = 0; i <= upper; ++i) {
        sum_s[i] += sum_l[i];
        if (i < upper) sum_l[i + 1] += sum_s[i];
    }
    std::vector<int> buf(upper + 1);
    auto induce = [&](const std::vector<int>& lms) {
        std::fill(sa.begin(), sa.end(), -1);
        std::copy(sum_s.begin(), sum_s.end(), buf.begin());
        for (int d : lms)
            if (d != n) sa[buf[s[d]]++] = d;
        std::copy(sum_l.begin(), sum_l.end(), buf.begin());
        sa[buf[s[n - 1]]++] = n - 1;
        for (int i = 0; i < n; ++i) {
            const int v = sa[i];
            if (v >= 1 && !ls[v - 1]) sa[buf[s[v - 1]]++] = v - 1;
        }
        std::copy(sum_l.begin(), sum_l.end(), buf.begin());
        for (int i = n - 1; i >= 0; --i) {
            const int v = sa[i];
            if (v >= 1 && ls[v - 1]) sa[--buf[s[v - 1] + 1]] = v - 1;
        }
    };
    std::vector<int> lms_map(n + 1, -1), lms;
    int m = 0;
    for (int i = 1; i < n; ++i)
        if (!ls[i - 1] && ls[i]) lms_map[i] = m++;
    lms.reserve(m);
    for (int i = 1; i < n; ++i)
        if (!ls[i - 1] && ls[i]) lms.push_back(i);
    induce(lms);
    if (m) {
        std::vector<int> sorted_lms;
        sorted_lms.reserve(m);
        for (int v : sa)
            if (lms_map[v] != -1) sorted_lms.push_back(v);
        std::vector<int> rec_s(m);
        int rec_upper = 0;
        rec_s[lms_map[sorted_lms[0]]] = 0;
        for (int i = 1; i < m; ++i) {
            int l = sorted_lms[i - 1], r = sorted_lms[i];
            const int end_l = (lms_map[l] + 1 < m) ? lms[lms_map[l] + 1] : n;
            const int end_r = (lms_map[r] + 1 < m) ? lms[lms_map[r] + 1] : n;
            bool same = true;
            if (end_l - l != end_r - r) {
                same = false;
            } else {
                while (l < end_l) {
                    if (s[l] != s[r]) break;
                    ++l, ++r;
                }
                if (l == n || s[l] != s[r]) same = false;
            }
            if (!same) ++rec_upper;
            rec_s[lms_map[sorted_lms[i]]] = rec_upper;
        }
        const std::vector<int> rec_sa = sa_is(rec_s, rec_upper);
        for (int i = 0; i < m; ++i) sorted_lms[i] = lms[rec_sa[i]];
        induce(sorted_lms);
    }
    return sa;
}
inline std::vector<int> lcp_kasai(const std::vector<int>& s, const std::vector<int>& sa) {
    const int n = (int)s.size();
    if (n == 0) return {};
    std::vector<int> rnk(n), lcp(n - 1);
    for (int i = 0; i < n; ++i) rnk[sa[i]] = i;
    int h = 0;
    for (int i = 0; i < n; ++i) {
        if (h > 0) --h;
        if (rnk[i] == 0) continue;
        const int j = sa[rnk[i] - 1];
        while (j + h < n && i + h < n && s[j + h] == s[i + h]) ++h;
        lcp[rnk[i] - 1] = h;
    }
    return lcp;
}
// 列 → 0〜upper の整数の列（char は 0〜255。数は圧縮）
template <class Seq>
std::vector<int> sa_ints(const Seq& s0, int& upper) {
    const auto& s = seq(s0);
    std::vector<int> r(s.size());
    if constexpr (std::is_same<elem_t<Seq>, char>::value || std::is_same<elem_t<Seq>, signed char>::value || std::is_same<elem_t<Seq>, unsigned char>::value) {
        for (std::size_t i = 0; i < s.size(); ++i) r[i] = (int)(unsigned char)s[i];
        upper = 255;
    } else {
        const CompressedValues c = compress_values(s);
        r = c.ids;
        upper = std::max(0, c.size() - 1);
    }
    return r;
}
}  // namespace str_internal

template <class Seq>
std::vector<int> make_suffix_array(const Seq& s) {  // 接尾辞を辞書順に並べた開始位置（SA-IS）  O(N)（数の列は圧縮で O(N log N)）
    int upper = 0;
    const auto v = str_internal::sa_ints(s, upper);
    return str_internal::sa_is(v, upper);
}
template <class Seq>
std::vector<int> make_lcp_array(const Seq& s, const std::vector<int>& sa) {  // lcp[i] = sa[i] と sa[i+1] から始まる接尾辞の LCP（長さ N - 1）  O(N)
    int upper = 0;
    return str_internal::lcp_kasai(str_internal::sa_ints(s, upper), sa);
}

// 2 か所の同じ部分文字列（長さ len と 2 つの開始位置。ないときは len = 0, pos = -1）
//   longest_repeated 系：同じ文字列の中の 2 か所（pos_a < pos_b）/ longest_common_substring：a での位置 pos_a と b での位置 pos_b
struct CommonSubstring {
    int len = 0;
    int pos_a = -1, pos_b = -1;
};

struct SuffixArray {
    std::vector<int> sa;       // sa[k] = 辞書順で k 番目の接尾辞の開始位置
    std::vector<int> rank;     // rank[p] = 位置 p から始まる接尾辞の順位（sa の逆）
    std::vector<int> lcp_adj;  // lcp_adj[k] = sa[k] と sa[k+1] の接尾辞の LCP（長さ N - 1）

    SuffixArray() = default;
    template <class Seq>
    explicit SuffixArray(const Seq& s0) {
        const auto& s = str_internal::seq(s0);
        n_ = (int)s.size();
        val_.resize(n_);
        for (int i = 0; i < n_; ++i) val_[i] = str_internal::code_of(s[i]);
        int upper = 0;
        const auto v = str_internal::sa_ints(s, upper);
        build(v, upper);
    }
    // 複数の文字列をまとめて（間に区切りを入れる。区切りはどの文字より小さく、互いに違う）。sa には区切りの位置（owner = -1）も入る
    template <class Seq>
    static SuffixArray of_strings(const std::vector<Seq>& vs) {
        SuffixArray R;
        const int K = (int)vs.size();
        long long total = K;
        for (const auto& x : vs) total += (long long)str_internal::seq(x).size();
        if (total > 2000000000LL) str_internal::fail("SuffixArray::of_strings", "長さの合計が大きすぎる");
        R.n_ = (int)total, R.multi_ = true;
        R.val_.resize(R.n_), R.lim_.resize(R.n_), R.owner_.resize(R.n_), R.start_.resize(K + 1);
        std::vector<long long> all;
        for (const auto& x : vs) {
            const auto& s = str_internal::seq(x);
            for (std::size_t i = 0; i < s.size(); ++i) all.push_back(str_internal::code_of(s[i]));
        }
        const CompressedValues cv = compress_values(all);
        std::vector<int> v(R.n_);
        int p = 0, q = 0;
        for (int t = 0; t < K; ++t) {
            const int len = (int)str_internal::seq(vs[t]).size();
            R.start_[t] = p;
            for (int i = 0; i < len; ++i, ++p, ++q) v[p] = K + cv.ids[q], R.val_[p] = all[q], R.lim_[p] = R.start_[t] + len, R.owner_[p] = t;
            v[p] = t, R.val_[p] = 0, R.lim_[p] = p, R.owner_[p] = -1, ++p;  // 区切り
        }
        R.start_[K] = p;
        R.build(v, K + std::max(0, cv.size() - 1));
        return R;
    }

    int size() const { return n_; }  // 全体の長さ（of_strings なら区切りも含む）
    int length_at(int p) const {  // 位置 p から始まる接尾辞の（自分の文字列の中での）長さ（0 ≤ p ≤ N。区切りと p = N は 0）
        check_pos("SuffixArray::length_at", p);
        return p == n_ ? 0 : lim(p) - p;
    }
    int owner(int p) const {  // 位置 p がどの文字列か（区切りは -1。0 ≤ p < N）
        check_index("SuffixArray::owner", p);
        return multi_ ? owner_[p] : 0;
    }
    int offset(int p) const {  // 位置 p がその文字列の何文字目か（区切りの位置では止める。owner(p) ≥ 0 を先に確かめる）
        check_index("SuffixArray::offset", p);
        if (!multi_) return p;
        if (owner_[p] < 0) str_internal::fail("SuffixArray::offset", "位置 " + std::to_string(p) + " は区切り（owner = -1）");
        return p - start_[owner_[p]];
    }
    int position(int t, int off) const {  // t 番目の文字列の off 文字目の位置（0 ≤ off ≤ その長さ。off = 長さ は直後の区切り）
        const int K = multi_ ? (int)start_.size() - 1 : 1;
        if (t < 0 || t >= K) str_internal::fail("SuffixArray::position", "文字列の番号 " + std::to_string(t) + " は 0〜" + std::to_string(K - 1) + " の外");
        const int len = multi_ ? start_[t + 1] - start_[t] - 1 : n_;
        if (off < 0 || off > len) str_internal::fail("SuffixArray::position", "off = " + std::to_string(off) + " は 0〜" + std::to_string(len) + " の外");
        return multi_ ? start_[t] + off : off;
    }
    long long value(int p) const {  // 位置 p の値（char は 0〜255。区切りは 0）
        check_index("SuffixArray::value", p);
        return val_[p];
    }

    // 位置 i, j から始まる接尾辞の LCP（0 ≤ i, j ≤ N）  O(1)（初回に疎テーブルを O(N log N) で作る）
    int lcp(int i, int j) const {
        check_pos("SuffixArray::lcp", i), check_pos("SuffixArray::lcp", j);
        if (i == n_ || j == n_) return 0;
        if (i == j) return length_at(i);
        int a = rank[i], b = rank[j];
        if (a > b) std::swap(a, b);
        return range_min(a, b);
    }
    // 部分文字列 s[l1, r1) と s[l2, r2) の辞書順（-1 / 0 / 1）  O(1)
    int compare(int l1, int r1, int l2, int r2) const {
        check_sub("SuffixArray::compare", l1, r1), check_sub("SuffixArray::compare", l2, r2);
        const int m1 = r1 - l1, m2 = r2 - l2, m = std::min(m1, m2);
        if (m == 0 || lcp(l1, l2) >= m) return m1 == m2 ? 0 : (m1 < m2 ? -1 : 1);
        return rank[l1] < rank[l2] ? -1 : 1;
    }
    // pat で始まる接尾辞の SA 上の区間 [l, r)（なければ空。pat が空なら全部 [0, N)）  O(M log N)
    template <class P>
    StrSpan find_range(const P& pat0) const {
        const auto& pat = str_internal::seq(pat0);
        std::vector<long long> p(pat.size());
        for (std::size_t i = 0; i < pat.size(); ++i) p[i] = str_internal::code_of(pat[i]);
        int lo = 0, hi = n_;  // 最初の「接尾辞 ≥ pat」
        while (lo < hi) {
            const int mid = (lo + hi) >> 1;
            if (cmp_pat(sa[mid], p) < 0) lo = mid + 1;
            else hi = mid;
        }
        int lo2 = lo, hi2 = n_;  // 最初の「接尾辞 > pat（pat で始まらない）」
        while (lo2 < hi2) {
            const int mid = (lo2 + hi2) >> 1;
            if (cmp_pat(sa[mid], p) <= 0) lo2 = mid + 1;
            else hi2 = mid;
        }
        return {lo, lo2};
    }
    template <class P>
    long long count(const P& pat) const {  // pat の出現回数（pat が空なら置ける位置の数：N + 1。of_strings では各文字列の長さ + 1 の和 = N）  O(M log N)
        if (str_internal::seq(pat).size() == 0) return empty_count();
        return find_range(pat).len();
    }
    template <class P>
    std::vector<int> find_all(const P& pat) const {  // pat の出現位置（昇順。pat が空なら 0〜N（of_strings では 0〜N-1））  O(M log N + 個数 log)
        if (str_internal::seq(pat).size() == 0) {
            std::vector<int> v((std::size_t)empty_count());
            for (std::size_t i = 0; i < v.size(); ++i) v[i] = (int)i;
            return v;
        }
        const StrSpan r = find_range(pat);
        std::vector<int> v(sa.begin() + r.l, sa.begin() + r.r);
        std::sort(v.begin(), v.end());
        return v;
    }
    long long count_occurrences(int l, int r) const {  // 部分文字列 s[l, r) の出現回数（自分も含む）  O(log N)
        check_sub("SuffixArray::count_occurrences", l, r);
        if (l == r) return empty_count();
        const StrSpan g = occurrence_range(l, r);
        return g.len();
    }
    StrSpan occurrence_range(int l, int r) const {  // s[l, r)（空でない）で始まる接尾辞の SA 上の区間  O(log N)
        check_sub("SuffixArray::occurrence_range", l, r);
        if (l == r) str_internal::fail("SuffixArray::occurrence_range", "空の部分文字列");
        const int m = r - l, k = rank[l];
        ensure_table();
        int lo = 0, hi = k;  // 最小の a で、a..k の LCP がすべて m 以上
        while (lo < hi) {
            const int mid = (lo + hi) >> 1;
            if (range_min(mid, k) >= m) hi = mid;
            else lo = mid + 1;
        }
        int lo2 = k, hi2 = n_ - 1;  // 最大の b で、k..b の LCP がすべて m 以上
        while (lo2 < hi2) {
            const int mid = (lo2 + hi2 + 1) >> 1;
            if (range_min(k, mid) >= m) lo2 = mid;
            else hi2 = mid - 1;
        }
        return {lo, lo2 + 1};
    }

    // ---- 異なる部分文字列 ----
    long long count_distinct_substrings() const {  // 異なる部分文字列（空でない）の数  O(1)（初回 O(N)）
        ensure_newcnt();
        return newcnt_.back();
    }
    std::vector<long long> distinct_substrings_by_length() const {  // r[len] = 長さ len の異なる部分文字列の数（len = 0〜最長）  O(N)
        int mx = 0;
        for (int p = 0; p < n_; ++p) mx = std::max(mx, length_at(p));
        std::vector<long long> d(mx + 2, 0);
        for (int k = 0; k < n_; ++k) {
            const int lo = (k ? lcp_adj[k - 1] : 0) + 1, hi = length_at(sa[k]);
            if (lo <= hi) ++d[lo], --d[hi + 1];
        }
        std::vector<long long> r(mx + 1, 0);
        for (int t = 1; t <= mx; ++t) d[t] += d[t - 1], r[t] = d[t];
        return r;
    }
    StrSpan kth_substring(long long k) const {  // 異なる部分文字列の辞書順 k 番目（0 から）  O(log N)
        ensure_newcnt();
        if (k < 0 || k >= newcnt_.back()) str_internal::fail("SuffixArray::kth_substring", "k は 0〜（異なる部分文字列の数 - 1）");
        const int idx = str_internal::upper_index(newcnt_, k) - 1;
        const int len = (idx ? lcp_adj[idx - 1] : 0) + (int)(k - newcnt_[idx]) + 1;
        return {sa[idx], sa[idx] + len};
    }
    long long substring_rank(int l, int r) const {  // 異なる部分文字列のうち s[l, r)（空でない）より小さいものの数  O(log N)
        const StrSpan g = occurrence_range(l, r);
        ensure_newcnt();
        const int k0 = g.l;
        return newcnt_[k0] + (r - l - 1 - (k0 ? lcp_adj[k0 - 1] : 0));
    }
    long long count_substrings_with_dup() const {  // 部分文字列（位置が違えば別）の数 = Σ 長さ  O(1)（初回 O(N)）
        ensure_blocks();
        return bsum_.back();
    }
    StrSpan kth_substring_with_dup(long long k) const {  // 位置が違えば別として数えたときの辞書順 k 番目（0 から）  O(log N)（初回 O(N)）
        ensure_blocks();
        if (k < 0 || k >= bsum_.back()) str_internal::fail("SuffixArray::kth_substring_with_dup", "k は 0〜（部分文字列の数 - 1）");
        const int b = str_internal::upper_index(bsum_, k) - 1;
        const Block& B = blocks_[b];
        const int len = B.t_lo + (int)((k - bsum_[b]) / B.mult);
        return {sa[B.k], sa[B.k] + len};
    }

    // ---- 繰り返し・LCP の和 ----
    CommonSubstring longest_repeated() const {  // 2 か所以上に現れる最長の部分文字列（重なってよい。同じ長さなら辞書順最小）  O(N)
        CommonSubstring c;
        for (int k = 0; k + 1 < n_; ++k)
            if (lcp_adj[k] > c.len) c.len = lcp_adj[k], c.pos_a = std::min(sa[k], sa[k + 1]), c.pos_b = std::max(sa[k], sa[k + 1]);
        return c;
    }
    CommonSubstring longest_repeated_non_overlapping() const {  // 重ならずに 2 か所に現れる最長の部分文字列  O(N log N)
        int mx = 0;
        for (int x : lcp_adj) mx = std::max(mx, x);
        auto test = [&](int L, CommonSubstring* out) {
            for (int k = 0; k < n_;) {  // lcp_adj ≥ L でつながる SA の区間ごとに、最小と最大の位置
                int j = k, mn = sa[k], mxp = sa[k];
                while (j + 1 < n_ && lcp_adj[j] >= L) ++j, mn = std::min(mn, sa[j]), mxp = std::max(mxp, sa[j]);
                if (j > k && mxp - mn >= L) {
                    if (out) *out = {L, mn, mxp};
                    return true;
                }
                k = j + 1;
            }
            return false;
        };
        int lo = 0, hi = mx;
        while (lo < hi) {
            const int mid = (lo + hi + 1) >> 1;
            if (test(mid, nullptr)) lo = mid;
            else hi = mid - 1;
        }
        CommonSubstring c;
        if (lo > 0) test(lo, &c);
        return c;
    }
    std::vector<long long> sum_lcp_with_all() const {  // r[i] = Σ_j lcp(i, j)（j はすべての位置。j = i も含む）  O(N)
        std::vector<long long> left(n_, 0), right(n_, 0), r(n_, 0);
        sweep(left, false), sweep(right, true);
        for (int k = 0; k < n_; ++k) r[sa[k]] = left[k] + right[k] + length_at(sa[k]);
        return r;
    }
    long long sum_lcp_all_pairs() const {  // Σ_{i < j} lcp(i, j)（すべての位置の組）  O(N)
        std::vector<long long> left(n_, 0);
        sweep(left, false);
        long long t = 0;
        for (long long x : left) t += x;
        return t;
    }

   private:
    struct Block {
        int k, t_lo, t_hi;
        long long mult;
    };
    int lim(int p) const { return multi_ ? lim_[p] : n_; }
    void check_pos(const char* func, int i) const {
        if (i < 0 || i > n_) str_internal::fail(func, "位置 " + std::to_string(i) + " は 0〜N（= " + std::to_string(n_) + "）の外");
    }
    void check_index(const char* func, int i) const {
        if (i < 0 || i >= n_) str_internal::fail(func, "位置 " + std::to_string(i) + " は 0〜N-1（N = " + std::to_string(n_) + "）の外");
    }
    long long empty_count() const { return multi_ ? n_ : (long long)n_ + 1; }  // 空の列を置ける位置の数
    void check_sub(const char* func, int l, int r) const {
        if (l < 0 || l > n_ || r < l || (l < n_ && r > lim(l)) || (l == n_ && r != l))
            str_internal::fail(func, "部分文字列 [" + std::to_string(l) + ", " + std::to_string(r) + ") が範囲外（of_strings では 1 つの文字列の中で）");
    }
    void build(const std::vector<int>& v, int upper) {
        sa = str_internal::sa_is(v, upper);
        rank.assign(n_, 0);
        for (int k = 0; k < n_; ++k) rank[sa[k]] = k;
        lcp_adj = str_internal::lcp_kasai(v, sa);
    }
    // 接尾辞 p と pat：p が pat より小さい（pat で始まらず小さい）なら -1、pat で始まるなら 0、大きいなら 1
    int cmp_pat(int p, const std::vector<long long>& pat) const {
        const int e = lim(p);
        for (std::size_t t = 0; t < pat.size(); ++t) {
            if (p + (int)t >= e) return -1;
            const long long x = val_[p + t];
            if (x != pat[t]) return x < pat[t] ? -1 : 1;
        }
        return 0;
    }
    void ensure_table() const {
        if (built_ || n_ <= 1) return;
        built_ = true;
        const int m = n_ - 1;
        int lg = 1;
        while ((1 << lg) <= m) ++lg;
        tab_.assign(lg, std::vector<int>());
        tab_[0] = lcp_adj;
        for (int j = 1; j < lg; ++j) {
            const int len = m - (1 << j) + 1;
            if (len <= 0) break;
            tab_[j].resize(len);
            for (int i = 0; i < len; ++i) tab_[j][i] = std::min(tab_[j - 1][i], tab_[j - 1][i + (1 << (j - 1))]);
        }
    }
    int range_min(int a, int b) const {  // min lcp_adj[a, b)（a < b）
        ensure_table();
        const int j = 31 - __builtin_clz((unsigned)(b - a));
        return std::min(tab_[j][a], tab_[j][b - (1 << j)]);
    }
    void ensure_newcnt() const {
        if (!newcnt_.empty()) return;
        newcnt_.assign(n_ + 1, 0);
        for (int k = 0; k < n_; ++k) newcnt_[k + 1] = newcnt_[k] + std::max(0, length_at(sa[k]) - (k ? lcp_adj[k - 1] : 0));
    }
    void ensure_blocks() const {
        if (!bsum_.empty()) return;
        // k を大きい方から。st は lcp_adj[k..] の「その位置までの最小」が下がる位置（値は狭義に減る）
        std::vector<std::pair<int, int>> st;  // (値, 位置)
        std::vector<Block> tmp;
        for (int k = n_ - 1; k >= 0; --k) {
            if (k <= n_ - 2) {
                while (!st.empty() && st.back().first >= lcp_adj[k]) st.pop_back();
                st.push_back({lcp_adj[k], k});
            }
            const int L = length_at(sa[k]), lp = k ? lcp_adj[k - 1] : 0;
            // t ∈ (v_0, L]：自分だけ（mult 1）。t ∈ (v_{i+1}, v_i]：mult = j_{i+1} - k + 1（なければ n - k）
            const int top = (int)st.size() - 1;
            if (k <= n_ - 2) {
                const int v0 = st[top].first;
                if (L > std::max(lp, v0)) tmp.push_back({k, std::max(lp, v0) + 1, L, 1});
                for (int i = top; i >= 0 && st[i].first > lp; --i) {
                    const int vi = st[i].first, vnext = i > 0 ? st[i - 1].first : 0;
                    const long long mult = (i > 0 ? st[i - 1].second : n_ - 1) - k + 1;
                    const int lo = std::max(lp, vnext) + 1;
                    if (lo <= vi) tmp.push_back({k, lo, vi, mult});
                }
            } else if (L > lp) {
                tmp.push_back({k, lp + 1, L, 1});
            }
        }
        // tmp は k の降順、同じ k の中では t の降順 → 反転すると k の昇順・t の昇順
        blocks_.assign(tmp.rbegin(), tmp.rend());
        bsum_.assign(blocks_.size() + 1, 0);
        for (std::size_t b = 0; b < blocks_.size(); ++b) bsum_[b + 1] = bsum_[b] + blocks_[b].mult * (blocks_[b].t_hi - blocks_[b].t_lo + 1);
    }
    void sweep(std::vector<long long>& out, bool from_right) const {  // out[k] = Σ_{k' < k}（右からなら k' > k）min lcp
        std::vector<std::pair<int, long long>> st;  // (値, 個数)
        long long cur = 0;
        for (int step = 1; step < n_; ++step) {
            const int k = from_right ? n_ - 1 - step : step;
            const int v = from_right ? lcp_adj[k] : lcp_adj[k - 1];
            long long c = 1;
            while (!st.empty() && st.back().first >= v) cur -= (long long)st.back().first * st.back().second, c += st.back().second, st.pop_back();
            st.push_back({v, c}), cur += (long long)v * c;
            out[k] = cur;
        }
    }
    int n_ = 0;
    bool multi_ = false;
    std::vector<long long> val_;
    std::vector<int> lim_, owner_, start_;
    mutable bool built_ = false;
    mutable std::vector<std::vector<int>> tab_;
    mutable std::vector<long long> newcnt_, bsum_;
    mutable std::vector<Block> blocks_;
};

// ---- 2 つ・K 個の文字列の最長共通部分文字列（連続する部分。飛び飛びの「部分列」は 06_lcs_edit） ----
template <class A, class B>
CommonSubstring longest_common_substring(const A& a, const B& b) {  // 長さと、a での位置 pos_a・b での位置 pos_b（なければ len = 0, -1）  O(N + M)
    using V = std::vector<long long>;
    const SuffixArray S = SuffixArray::of_strings(std::vector<V>{str_internal::codes_of(a), str_internal::codes_of(b)});
    CommonSubstring c;
    for (int k = 0; k + 1 < S.size(); ++k) {
        const int p = S.sa[k], q = S.sa[k + 1];
        if (S.owner(p) < 0 || S.owner(q) < 0 || S.owner(p) == S.owner(q)) continue;
        if (S.lcp_adj[k] > c.len) {
            c.len = S.lcp_adj[k];
            c.pos_a = S.offset(S.owner(p) == 0 ? p : q), c.pos_b = S.offset(S.owner(p) == 0 ? q : p);
        }
    }
    return c;
}
template <class A, class B>
std::vector<StrSpan> longest_common_substrings_all(const A& a, const B& b) {  // 最長共通部分文字列をすべて（異なるもの・辞書順。a での位置）  O(N + M)
    using V = std::vector<long long>;
    const SuffixArray S = SuffixArray::of_strings(std::vector<V>{str_internal::codes_of(a), str_internal::codes_of(b)});
    int best = 0;
    for (int k = 0; k + 1 < S.size(); ++k) {
        const int p = S.sa[k], q = S.sa[k + 1];
        if (S.owner(p) >= 0 && S.owner(q) >= 0 && S.owner(p) != S.owner(q)) best = std::max(best, S.lcp_adj[k]);
    }
    std::vector<StrSpan> r;
    if (best == 0) return r;
    int last = -1;  // 直前に入れたものの位置（全体の位置）
    for (int k = 0; k + 1 < S.size(); ++k) {
        const int p = S.sa[k], q = S.sa[k + 1];
        if (S.owner(p) < 0 || S.owner(q) < 0 || S.owner(p) == S.owner(q) || S.lcp_adj[k] != best) continue;
        const int pa = S.owner(p) == 0 ? p : q;
        if (last >= 0 && S.lcp(last, pa) >= best) continue;
        last = pa;
        r.push_back({S.offset(pa), S.offset(pa) + best});
    }
    return r;
}
template <class Seq>
StrSpan longest_common_substring_k(const std::vector<Seq>& vs) {  // すべての文字列に現れる最長の部分文字列（vs[0] での位置。なければ {0, 0}）  O(全体の長さ)
    const int K = (int)vs.size();
    if (K == 0) return {0, 0};
    if (K == 1) return {0, (int)str_internal::seq(vs[0]).size()};
    const SuffixArray S = SuffixArray::of_strings(vs);
    const int n = S.size();
    std::vector<int> cnt(K, 0);
    int kinds = 0, best = 0, best_pos = -1, last0 = -1;  // last0 = 窓の中で最後の「0 番目の文字列」の SA の位置
    std::vector<int> dq(n + 1);  // lcp_adj の位置の単調な列（最小を先頭に）
    int qh = 0, qt = 0;
    for (int l = 0, r = 0; r < n; ++r) {
        const int o = S.owner(S.sa[r]);
        if (o >= 0 && cnt[o]++ == 0) ++kinds;
        if (o == 0) last0 = r;
        if (r > 0) {
            while (qt > qh && S.lcp_adj[dq[qt - 1]] >= S.lcp_adj[r - 1]) --qt;
            dq[qt++] = r - 1;
        }
        while (kinds == K) {
            if (r > l) {
                while (qh < qt && dq[qh] < l) ++qh;
                const int v = S.lcp_adj[dq[qh]];
                if (v > best) best = v, best_pos = S.offset(S.sa[last0]);
            }
            const int o2 = S.owner(S.sa[l]);
            if (o2 >= 0 && --cnt[o2] == 0) --kinds;
            ++l;
        }
    }
    if (best == 0) return {0, 0};
    return {best_pos, best_pos + best};
}
#endif  // GG_SUFFIX_ARRAY

#ifndef GG_STR_BASE
#define GG_STR_BASE
// ---- 文字列基盤（文字列第二世代）  文字の種類 Alphabet・文字の数と位置・ランレングス・変換・比較・区切り ----
// 列（string / vector<整数>）を受け取る関数はどちらでもよい。文字列リテラル（"abc"）も渡せる
// 範囲は半開区間 [l, r)（0 始まり）。位置が「ない」ときは -1。使える範囲の外の引数はメッセージを出して止める
namespace str_internal {
[[noreturn]] inline void fail(const char* func, const std::string& msg) {
    std::cerr << func << ": " << msg << "\n";
    std::abort();
}
// 文字列リテラルを std::string に（それ以外はそのまま参照を返す）
template <class T>
const T& seq(const T& s) { return s; }
template <std::size_t N>
std::string seq(const char (&s)[N]) { return std::string(s); }
inline std::string seq(const char* s) { return std::string(s); }
template <class T>
using seq_t = std::decay_t<decltype(seq(std::declval<const T&>()))>;
template <class T>
using elem_t = std::decay_t<decltype(std::declval<const seq_t<T>&>()[0])>;
// 要素を「値」に（char は unsigned char として 0〜255。負の添字にしない）
template <class T>
long long code_of(T x) {
    if constexpr (std::is_same<T, char>::value || std::is_same<T, signed char>::value || std::is_same<T, unsigned char>::value)
        return (long long)(unsigned char)x;
    else
        return (long long)x;
}
// 範囲の確認（止めるときの文字列作りは別の関数に分けて、確認そのものを軽くする）
[[noreturn]] __attribute__((noinline, cold)) inline void range_fail(const char* func, long long l, long long r, long long n) {
    fail(func, "範囲 [l, r) は 0 ≤ l ≤ r ≤ 長さ で指定する（l = " + std::to_string(l) + ", r = " + std::to_string(r) + ", 長さ = " + std::to_string(n) + "）");
}
inline void check_range(const char* func, long long l, long long r, long long n) {
    if (__builtin_expect(l < 0 || r > n || l > r, 0)) range_fail(func, l, r, n);
}
// a[lo, hi) の中で x 以上の最初の位置（a は昇順）。分岐の少ない二分探索
inline int lower_index(const int* a, int lo, int hi, int x) {
    int len = hi - lo;
    if (len <= 0) return lo;
    const int* b = a + lo;
    while (len > 1) {
        const int half = len >> 1;
        b += (b[half - 1] < x) ? half : 0;
        len -= half;
    }
    return (int)(b - a) + (*b < x);
}
// 昇順の a の中で x より大きい最初の位置（std::upper_bound と同じ。デバッグビルドでも O(log N)）
template <class T>
int upper_index(const std::vector<T>& a, const T& x) {
    int lo = 0, hi = (int)a.size();
    while (lo < hi) {
        const int mid = (int)((unsigned)(lo + hi) >> 1);
        if (!(x < a[mid])) lo = mid + 1;
        else hi = mid;
    }
    return lo;
}
// 数え上げの足し算：T が整数型なら溢れたら止める（mint・double などはそのまま）
template <class T>
T checked_add(const T& a, const T& b, const char* func) {
    if constexpr (std::is_integral<T>::value) {
        T r;
        if (__builtin_add_overflow(a, b, &r)) fail(func, "値が型に収まらない（答えを mod で求めるなら math_v2 の ModInt998 などを使う）");
        return r;
    } else {
        return a + b;
    }
}
// 数え上げの途中の足し算：T が整数型なら溢れたことを印 oa に残すだけ（途中で溢れても答えに使わなければ止めない。最後に印を見て止める）
template <class T>
void add_count(T& a, bool& oa, const T& b, bool ob) {
    if constexpr (std::is_integral<T>::value) {
        T r;
        oa = __builtin_add_overflow(a, b, &r) || oa || ob;
        a = r;
    } else {
        a = a + b;
    }
}
[[noreturn]] inline void count_overflow(const char* func) { fail(func, "値が型に収まらない（答えを mod で求めるなら math_v2 の ModInt998 などを使う）"); }
// 列 → 値の列（char は 0〜255）
template <class Seq>
std::vector<long long> codes_of(const Seq& s0) {
    const auto& s = seq(s0);
    std::vector<long long> c(s.size());
    for (std::size_t i = 0; i < s.size(); ++i) c[i] = code_of(s[i]);
    return c;
}
// 辞書順の比較（-1 / 0 / 1）。char は unsigned char として比べる（std::string と同じ順）
template <class A, class B>
int lex_compare(const A& a, const B& b) {
    const std::size_t n = a.size(), m = b.size(), k = n < m ? n : m;
    for (std::size_t i = 0; i < k; ++i) {
        const long long x = code_of(a[i]), y = code_of(b[i]);
        if (x != y) return x < y ? -1 : 1;
    }
    return n == m ? 0 : (n < m ? -1 : 1);
}
}  // namespace str_internal

// ---- 部分文字列の位置 [l, r) ----
struct StrSpan {
    int l = 0, r = 0;
    int len() const { return r - l; }
    bool empty() const { return r <= l; }
    friend bool operator==(const StrSpan& a, const StrSpan& b) { return a.l == b.l && a.r == b.r; }
    friend bool operator!=(const StrSpan& a, const StrSpan& b) { return !(a == b); }
    friend std::ostream& operator<<(std::ostream& os, const StrSpan& sp) { return os << "[" << sp.l << ", " << sp.r << ")"; }  // デバッグ用
};
template <class Seq>
str_internal::seq_t<Seq> substr_of(const Seq& s0, StrSpan sp) {  // s[sp.l, sp.r) のコピー（string なら string、vector なら vector）
    const auto& s = str_internal::seq(s0);
    str_internal::check_range("substr_of", sp.l, sp.r, (long long)s.size());
    return str_internal::seq_t<Seq>(s.begin() + sp.l, s.begin() + sp.r);
}

// ---- 文字の種類（Alphabet）：文字・値 → 番号 0〜size()-1。番号の順が辞書順として扱われる ----
//   Alphabet::lower()（a〜z。既定）/ upper()（A〜Z）/ digits()（0〜9）/ letters()（A〜Z, a〜z の 52 種。ASCII 順）
//   alnum()（0〜9, A〜Z, a〜z の 62 種。ASCII 順）/ ascii()（0〜127）/ bytes()（0〜255。日本語など 128 以上の char も）
//   of("()[]")（並べた文字だけ。並べた順が番号の順 = 辞書順）/ ints(K)（整数 0〜K-1）/ ints(lo, hi)（整数 [lo, hi)）
//   範囲外の文字・値を id に渡すとメッセージを出して止まる（contains で先に確かめられる）
struct Alphabet {
    Alphabet() : Alphabet(of("abcdefghijklmnopqrstuvwxyz")) {}
    static Alphabet lower() { return of("abcdefghijklmnopqrstuvwxyz"); }
    static Alphabet upper() { return of("ABCDEFGHIJKLMNOPQRSTUVWXYZ"); }
    static Alphabet digits() { return of("0123456789"); }
    static Alphabet letters() { return of("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"); }
    static Alphabet alnum() { return of("0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"); }
    static Alphabet ascii() { return code_range(0, 128); }
    static Alphabet bytes() { return code_range(0, 256); }
    static Alphabet of(const std::string& cs) {  // 並べた文字（同じ文字を 2 回並べると止まる）
        if (cs.empty()) str_internal::fail("Alphabet::of", "文字が 1 つもない");
        Alphabet A(0);
        A.kind_ = CHARS;
        for (std::size_t i = 0; i < cs.size(); ++i) {
            const unsigned char c = (unsigned char)cs[i];
            if (A.table_[c] >= 0) str_internal::fail("Alphabet::of", std::string("同じ文字が 2 回ある: '") + cs[i] + "'");
            A.table_[c] = (short)i;
            A.chars_[i] = c;
        }
        A.size_ = (int)cs.size();
        return A;
    }
    static Alphabet ints(long long K) { return ints(0, K); }  // 整数 0〜K-1
    static Alphabet ints(long long lo, long long hi) {         // 整数 [lo, hi)
        if (hi <= lo || (unsigned long long)hi - (unsigned long long)lo > 1000000000ULL)
            str_internal::fail("Alphabet::ints", "lo < hi かつ hi - lo ≤ 10^9 で指定する（大きな値は compress_values で 0〜K-1 に）");
        Alphabet A(0);
        A.kind_ = INTS;
        A.lo_ = lo;
        A.size_ = (int)(hi - lo);
        return A;
    }

    int size() const { return size_; }
    bool is_chars() const { return kind_ == CHARS; }  // 文字の表で決めたものか（of・lower など）。false なら ints
    template <class T>
    bool contains(T x) const { return index_of(str_internal::code_of(x)) >= 0; }
    template <class T>
    int id(T x) const {  // 文字・値の番号（範囲外は止める）
        const long long v = str_internal::code_of(x);
        const int r = index_of(v);
        if (r < 0) out_of_range(v);
        return r;
    }
    long long value(int i) const {  // 番号 i の値（文字なら文字コード 0〜255）
        check_index(i);
        return kind_ == CHARS ? (long long)chars_[i] : lo_ + i;
    }
    char ch(int i) const {  // 番号 i の文字（ints で 0〜255 の外の値なら止める）
        const long long v = value(i);
        if (v < 0 || v > 255) str_internal::fail("Alphabet::ch", "番号 " + std::to_string(i) + " の値 " + std::to_string(v) + " は文字にできない");
        return (char)(unsigned char)v;
    }
    std::string chars() const {  // 文字を番号の順に並べたもの（ints なら 0〜255 の範囲のときだけ）
        std::string s;
        for (int i = 0; i < size_; ++i) s += ch(i);
        return s;
    }

   private:
    enum Kind : unsigned char { CHARS, INTS };
    explicit Alphabet(int) { table_.fill(-1); }
    static Alphabet code_range(int lo, int hi) {
        Alphabet A(0);
        A.kind_ = CHARS;
        for (int c = lo; c < hi; ++c) A.table_[c] = (short)(c - lo), A.chars_[c - lo] = (unsigned char)c;
        A.size_ = hi - lo;
        return A;
    }
    int index_of(long long v) const {
        if (kind_ == CHARS) return (unsigned long long)v < 256ULL ? table_[(std::size_t)v] : -1;
        const unsigned long long d = (unsigned long long)v - (unsigned long long)lo_;
        return d < (unsigned long long)size_ ? (int)d : -1;
    }
    void check_index(int i) const {
        if (i < 0 || i >= size_) str_internal::fail("Alphabet", "番号 " + std::to_string(i) + " は 0〜" + std::to_string(size_ - 1) + " の外");
    }
    [[noreturn]] void out_of_range(long long v) const {
        std::string what = "値 " + std::to_string(v) + " は";
        if (kind_ == CHARS && v >= 32 && v < 127) what = std::string("文字 '") + (char)v + "'（コード " + std::to_string(v) + "）は";
        str_internal::fail("Alphabet::id", what + "文字の種類に入っていない（種類の数 " + std::to_string(size_) + "）");
    }
    Kind kind_ = CHARS;
    int size_ = 0;
    long long lo_ = 0;
    std::array<short, 256> table_{};
    std::array<unsigned char, 256> chars_{};
};

template <class Seq>
std::vector<int> to_ids(const Seq& s0, const Alphabet& A = Alphabet::lower()) {  // 各文字の番号の列
    const auto& s = str_internal::seq(s0);
    std::vector<int> r(s.size());
    for (std::size_t i = 0; i < s.size(); ++i) r[i] = A.id(s[i]);
    return r;
}
inline std::string from_ids(const std::vector<int>& ids, const Alphabet& A = Alphabet::lower()) {  // 番号の列 → 文字列
    std::string s(ids.size(), '\0');
    for (std::size_t i = 0; i < ids.size(); ++i) s[i] = A.ch(ids[i]);
    return s;
}

// ---- 数の列を 0〜K-1 に圧縮（大小の順を保つ）：Alphabet::ints(K) と組み合わせる ----
struct CompressedValues {
    std::vector<int> ids;            // 各要素の番号（0〜K-1。値の小さい順）
    std::vector<long long> values;   // 番号 → 元の値（昇順。K = values.size()）
    int size() const { return (int)values.size(); }
    Alphabet alphabet() const { return Alphabet::ints(values.empty() ? 1 : (long long)values.size()); }
};
template <class Seq>
CompressedValues compress_values(const Seq& s0) {
    const auto& s = str_internal::seq(s0);
    const std::size_t n = s.size();
    std::vector<std::pair<long long, int>> vi(n);
    for (std::size_t i = 0; i < n; ++i) vi[i] = {str_internal::code_of(s[i]), (int)i};
    std::sort(vi.begin(), vi.end());
    CompressedValues c;
    c.ids.resize(n);
    for (std::size_t i = 0; i < n; ++i) {
        if (i == 0 || vi[i].first != vi[i - 1].first) c.values.push_back(vi[i].first);
        c.ids[vi[i].second] = (int)c.values.size() - 1;
    }
    return c;
}

// ---- 文字の数 ----
template <class Seq>
std::vector<int> count_chars(const Seq& s0, const Alphabet& A = Alphabet::lower()) {  // 番号ごとの個数（長さ A.size()）
    const auto& s = str_internal::seq(s0);
    std::vector<int> c(A.size(), 0);
    for (std::size_t i = 0; i < s.size(); ++i) ++c[A.id(s[i])];
    return c;
}

// ---- 区間の文字の数・次 / 前の出現を O(1) で（前処理 O(Nσ)・メモリ (N+1)σ 個の int。N = 2×10^5, σ = 26 で約 21MB） ----
//   メモリを抑えたいとき・σ が大きいときは CharPositions（O(log N)・メモリ O(N + σ)）
struct CharPrefix {
    CharPrefix() = default;
    template <class Seq>
    explicit CharPrefix(const Seq& s0, const Alphabet& A = Alphabet::lower()) : A_(A) {
        const auto& s = str_internal::seq(s0);
        n_ = (int)s.size(), sg_ = A_.size();
        if ((long long)(n_ + 1) * sg_ > 600000000LL) str_internal::fail("CharPrefix", "(N + 1)σ が 6×10^8 を超える（メモリが足りない。CharPositions を使う）");
        cnt_.assign((std::size_t)(n_ + 1) * sg_, 0);
        std::vector<int> id(n_);
        start_.assign(sg_ + 1, 0);
        for (int i = 0; i < n_; ++i) {
            const int* src = cnt_.data() + (std::size_t)i * sg_;
            int* dst = cnt_.data() + (std::size_t)(i + 1) * sg_;
            std::copy(src, src + sg_, dst);
            id[i] = A_.id(s[i]);
            ++dst[id[i]], ++start_[id[i] + 1];
        }
        for (int k = 0; k < sg_; ++k) start_[k + 1] += start_[k];
        pos_.resize(n_);
        std::vector<int> cur(start_.begin(), start_.end() - 1);
        for (int i = 0; i < n_; ++i) pos_[cur[id[i]]++] = i;
    }
    int size() const { return n_; }
    const Alphabet& alphabet() const { return A_; }
    template <class C>
    int count(int l, int r, C c) const {  // [l, r) にある c の数  O(1)
        str_internal::check_range("CharPrefix::count", l, r, n_);
        const int k = A_.id(c);
        return cnt_[(std::size_t)r * sg_ + k] - cnt_[(std::size_t)l * sg_ + k];
    }
    template <class C>
    bool has(int l, int r, C c) const { return count(l, r, c) > 0; }  // [l, r) に c があるか  O(1)
    std::vector<int> counts(int l, int r) const {  // [l, r) の番号ごとの数  O(σ)
        str_internal::check_range("CharPrefix::counts", l, r, n_);
        std::vector<int> c(sg_);
        for (int k = 0; k < sg_; ++k) c[k] = cnt_[(std::size_t)r * sg_ + k] - cnt_[(std::size_t)l * sg_ + k];
        return c;
    }
    int distinct(int l, int r) const {  // [l, r) にある文字の種類の数  O(σ)
        str_internal::check_range("CharPrefix::distinct", l, r, n_);
        int d = 0;
        for (int k = 0; k < sg_; ++k) d += cnt_[(std::size_t)r * sg_ + k] != cnt_[(std::size_t)l * sg_ + k];
        return d;
    }
    template <class C>
    int next(int i, C c) const {  // i 以上で最初の c の位置（なければ -1。0 ≤ i ≤ N）  O(1)
        check_pos("CharPrefix::next", i);
        const int k = A_.id(c), m = cnt_[(std::size_t)i * sg_ + k];
        return start_[k] + m < start_[k + 1] ? pos_[start_[k] + m] : -1;
    }
    template <class C>
    int prev(int i, C c) const {  // i より前（[0, i)）で最後の c の位置（なければ -1。0 ≤ i ≤ N）  O(1)
        check_pos("CharPrefix::prev", i);
        const int k = A_.id(c), m = cnt_[(std::size_t)i * sg_ + k];
        return m > 0 ? pos_[start_[k] + m - 1] : -1;
    }
    template <class C>
    int kth(C c, long long k) const {  // c の k 番目（0 から）の位置（なければ -1）  O(1)
        const int id = A_.id(c);
        return (k >= 0 && k < start_[id + 1] - start_[id]) ? pos_[start_[id] + k] : -1;
    }
    template <class C>
    int total(C c) const {  // 全体での c の数  O(1)
        const int k = A_.id(c);
        return start_[k + 1] - start_[k];
    }

   private:
    void check_pos(const char* func, int i) const {
        if (i < 0 || i > n_) str_internal::fail(func, "位置 " + std::to_string(i) + " は 0〜N（= " + std::to_string(n_) + "）の外");
    }
    Alphabet A_;
    int n_ = 0, sg_ = 0;
    std::vector<int> cnt_;           // cnt_[i * σ + k] = [0, i) にある番号 k の数
    std::vector<int> start_, pos_;   // 番号ごとの出現位置（pos_[start_[k] .. start_[k+1])）
};

// ---- 文字の出現位置（前処理 O(N + σ)・メモリ O(N + σ)。σ が大きくても使える） ----
struct CharPositions {
    CharPositions() = default;
    template <class Seq>
    explicit CharPositions(const Seq& s0, const Alphabet& A = Alphabet::lower()) : A_(A) {
        const auto& s = str_internal::seq(s0);
        n_ = (int)s.size();
        const int sg = A_.size();
        std::vector<int> id(n_);
        start_.assign(sg + 1, 0);
        for (int i = 0; i < n_; ++i) id[i] = A_.id(s[i]), ++start_[id[i] + 1];
        for (int k = 0; k < sg; ++k) start_[k + 1] += start_[k];
        pos_.resize(n_), occ_.resize(n_);
        std::vector<int> cur(start_.begin(), start_.end() - 1);
        for (int i = 0; i < n_; ++i) occ_[i] = cur[id[i]] - start_[id[i]], pos_[cur[id[i]]++] = i;
    }
    int size() const { return n_; }
    const Alphabet& alphabet() const { return A_; }
    template <class C>
    int next(int i, C c) const {  // i 以上で最初の c の位置（なければ -1。0 ≤ i ≤ N）  O(log N)
        check_pos("CharPositions::next", i);
        const int k = A_.id(c), j = str_internal::lower_index(pos_.data(), start_[k], start_[k + 1], i);
        return j < start_[k + 1] ? pos_[j] : -1;
    }
    template <class C>
    int prev(int i, C c) const {  // i より前（[0, i)）で最後の c の位置（なければ -1。0 ≤ i ≤ N）  O(log N)
        check_pos("CharPositions::prev", i);
        const int k = A_.id(c), j = str_internal::lower_index(pos_.data(), start_[k], start_[k + 1], i);
        return j > start_[k] ? pos_[j - 1] : -1;
    }
    template <class C>
    int kth(C c, long long k) const {  // c の k 番目（0 から）の位置（なければ -1）  O(1)
        const int id = A_.id(c);
        return (k >= 0 && k < start_[id + 1] - start_[id]) ? pos_[start_[id] + k] : -1;
    }
    template <class C>
    int count(int l, int r, C c) const {  // [l, r) にある c の数  O(log N)
        str_internal::check_range("CharPositions::count", l, r, n_);
        const int k = A_.id(c);
        return str_internal::lower_index(pos_.data(), start_[k], start_[k + 1], r) - str_internal::lower_index(pos_.data(), start_[k], start_[k + 1], l);
    }
    template <class C>
    bool has(int l, int r, C c) const { return count(l, r, c) > 0; }
    template <class C>
    int total(C c) const {  // 全体での c の数  O(1)
        const int k = A_.id(c);
        return start_[k + 1] - start_[k];
    }
    template <class C>
    std::vector<int> positions(C c) const {  // c の位置の一覧（昇順）
        const int k = A_.id(c);
        return std::vector<int>(pos_.begin() + start_[k], pos_.begin() + start_[k + 1]);
    }
    int occurrence_index(int i) const {  // s[i] が「その文字の何番目（0 から）の出現」か  O(1)
        if (i < 0 || i >= n_) str_internal::fail("CharPositions::occurrence_index", "位置 " + std::to_string(i) + " は 0〜N-1 の外");
        return occ_[i];
    }

   private:
    void check_pos(const char* func, int i) const {
        if (i < 0 || i > n_) str_internal::fail(func, "位置 " + std::to_string(i) + " は 0〜N（= " + std::to_string(n_) + "）の外");
    }
    Alphabet A_;
    int n_ = 0;
    std::vector<int> start_, pos_, occ_;
};

// ---- ランレングス圧縮：連続する同じ値を {value, len, start} にまとめる ----
template <class T>
struct RunBlock {
    T value{};
    long long len = 0;    // 続く長さ
    long long start = 0;  // 元の列での開始位置
    friend std::ostream& operator<<(std::ostream& os, const RunBlock& b) { return os << "{" << b.value << ", " << b.len << ", " << b.start << "}"; }  // デバッグ用
};
template <class Seq>
std::vector<RunBlock<str_internal::elem_t<Seq>>> run_length(const Seq& s0) {
    const auto& s = str_internal::seq(s0);
    std::vector<RunBlock<str_internal::elem_t<Seq>>> r;
    std::size_t cnt = s.empty() ? 0 : 1;
    for (std::size_t i = 1; i < s.size(); ++i) cnt += !(s[i] == s[i - 1]);
    r.reserve(cnt);
    for (std::size_t i = 0, n = s.size(); i < n;) {
        std::size_t j = i + 1;
        while (j < n && s[j] == s[i]) ++j;
        r.push_back({s[i], (long long)(j - i), (long long)i});
        i = j;
    }
    return r;
}
template <class T>
auto run_length_decode(const std::vector<RunBlock<T>>& runs) {  // 元に戻す（T = char なら string、それ以外は vector<T>）
    using R = std::conditional_t<std::is_same<T, char>::value, std::string, std::vector<T>>;
    long long total = 0;
    for (const auto& b : runs) {
        if (b.len < 0) str_internal::fail("run_length_decode", "len が負");
        total += b.len;
        if (total > 2000000000LL) str_internal::fail("run_length_decode", "元の長さが 2×10^9 を超える");
    }
    R r;
    r.reserve((std::size_t)total);
    for (const auto& b : runs) r.insert(r.end(), (std::size_t)b.len, b.value);
    return r;
}

// ---- 変換 ----
inline std::string to_lower(std::string s) {  // 英大文字を小文字に（ほかはそのまま）
    for (char& c : s)
        if ('A' <= c && c <= 'Z') c = (char)(c - 'A' + 'a');
    return s;
}
inline std::string to_upper(std::string s) {
    for (char& c : s)
        if ('a' <= c && c <= 'z') c = (char)(c - 'a' + 'A');
    return s;
}
inline std::string swap_case(std::string s) {  // 大文字と小文字を入れ替える
    for (char& c : s) {
        if ('A' <= c && c <= 'Z') c = (char)(c - 'A' + 'a');
        else if ('a' <= c && c <= 'z') c = (char)(c - 'a' + 'A');
    }
    return s;
}
template <class Seq>
str_internal::seq_t<Seq> caesar_shift(const Seq& s0, long long k, const Alphabet& A = Alphabet::lower()) {  // A の文字を番号で k ずらす（負も。A にない文字はそのまま）
    auto s = str_internal::seq_t<Seq>(str_internal::seq(s0));
    const long long sg = A.size(), sh = ((k % sg) + sg) % sg;
    for (auto& c : s)
        if (A.contains(c)) c = (std::decay_t<decltype(c)>)A.value((int)((A.id(c) + sh) % sg));
    return s;
}
// 文字の置き換えを重ねる：replace(c, d) で「今 c になっている文字を全部 d に」。apply(s) で s に当てはめる（ABC342 C）
struct CharMap {
    explicit CharMap(const Alphabet& A = Alphabet::lower()) : A_(A), f_(A.size()) {
        for (int i = 0; i < A_.size(); ++i) f_[i] = i;
    }
    template <class C, class D>
    void replace(C c, D d) {  // 今 c になっているものを全部 d に  O(σ)
        const int x = A_.id(c), y = A_.id(d);
        if (x == y) return;
        for (int& v : f_)
            if (v == x) v = y;
    }
    template <class C, class D>
    void set(C c, D d) { f_[A_.id(c)] = A_.id(d); }  // 元の文字 c の行き先を d にする（ほかは変えない）
    template <class C>
    auto operator()(C c) const {  // 元の文字 c の行き先（A にない文字はそのまま）
        return A_.contains(c) ? (C)A_.value(f_[A_.id(c)]) : c;
    }
    template <class Seq>
    str_internal::seq_t<Seq> apply(const Seq& s0) const {  // s の各文字を行き先に（A にない文字はそのまま）
        auto s = str_internal::seq_t<Seq>(str_internal::seq(s0));
        for (auto& c : s) c = (*this)(c);
        return s;
    }
    CharMap compose(const CharMap& g) const {  // 「この置き換えのあとに g」を 1 つにしたもの
        CharMap h(A_);
        for (int i = 0; i < A_.size(); ++i) h.f_[i] = g.f_[f_[i]];
        return h;
    }
    const std::vector<int>& table() const { return f_; }  // 番号 → 行き先の番号

   private:
    Alphabet A_;
    std::vector<int> f_;
};
template <class Seq>
str_internal::seq_t<Seq> reversed(const Seq& s0) {  // 反転したコピー
    const auto& s = str_internal::seq(s0);
    return str_internal::seq_t<Seq>(s.rbegin(), s.rend());
}
template <class Seq>
str_internal::seq_t<Seq> min_with_reverse(const Seq& s0) {  // s と反転のうち辞書順で小さい方（向きを区別しないときの代表。ABC310 C）
    auto s = str_internal::seq_t<Seq>(str_internal::seq(s0));
    auto t = reversed(s);
    return str_internal::lex_compare(t, s) < 0 ? t : s;
}
template <class Seq>
str_internal::seq_t<Seq> sorted_chars(const Seq& s0) {  // 並べ替えたもの（並べ替えて同じかの比較用）
    auto s = str_internal::seq_t<Seq>(str_internal::seq(s0));
    std::sort(s.begin(), s.end(), [](const auto& a, const auto& b) { return str_internal::code_of(a) < str_internal::code_of(b); });
    return s;
}
template <class Seq>
str_internal::seq_t<Seq> repeat_string(const Seq& s0, long long k) {  // s を k 回つなげる
    const auto& s = str_internal::seq(s0);
    if (k < 0) str_internal::fail("repeat_string", "k ≥ 0 で使う");
    if (!s.empty() && k > 2000000000LL / (long long)s.size()) str_internal::fail("repeat_string", "長さが 2×10^9 を超える");
    str_internal::seq_t<Seq> r;
    r.reserve((std::size_t)(k * (long long)s.size()));
    for (long long i = 0; i < k; ++i) r.insert(r.end(), s.begin(), s.end());
    return r;
}

// ---- 比較 ----
template <class A, class B>
int common_prefix_length(const A& a0, const B& b0) {  // 共通の接頭辞の長さ
    const auto& a = str_internal::seq(a0);
    const auto& b = str_internal::seq(b0);
    const std::size_t k = a.size() < b.size() ? a.size() : b.size();
    std::size_t i = 0;
    while (i < k && str_internal::code_of(a[i]) == str_internal::code_of(b[i])) ++i;
    return (int)i;
}
template <class A, class B>
int common_suffix_length(const A& a0, const B& b0) {  // 共通の接尾辞の長さ
    const auto& a = str_internal::seq(a0);
    const auto& b = str_internal::seq(b0);
    const std::size_t k = a.size() < b.size() ? a.size() : b.size();
    std::size_t i = 0;
    while (i < k && str_internal::code_of(a[a.size() - 1 - i]) == str_internal::code_of(b[b.size() - 1 - i])) ++i;
    return (int)i;
}
template <class A, class B>
int hamming_distance(const A& a0, const B& b0) {  // 同じ位置で違う文字の数（長さが同じこと）
    const auto& a = str_internal::seq(a0);
    const auto& b = str_internal::seq(b0);
    if (a.size() != b.size()) str_internal::fail("hamming_distance", "長さが違う");
    int d = 0;
    for (std::size_t i = 0; i < a.size(); ++i) d += str_internal::code_of(a[i]) != str_internal::code_of(b[i]);
    return d;
}
template <class T, class S>
bool is_subsequence(const T& t0, const S& s0) {  // t が s の部分列（飛び飛びでよい）か  O(|s|)。何度も聞くなら 06_subsequence の NextTable
    const auto& t = str_internal::seq(t0);
    const auto& s = str_internal::seq(s0);
    std::size_t j = 0;
    for (std::size_t i = 0; i < s.size() && j < t.size(); ++i)
        if (str_internal::code_of(s[i]) == str_internal::code_of(t[j])) ++j;
    return j == t.size();
}
template <class Seq>
bool is_palindrome(const Seq& s0) {  // 全体が回文か  O(N)。区間を何度も聞くなら 04_palindrome の Palindromes
    const auto& s = str_internal::seq(s0);
    for (std::size_t i = 0, n = s.size(); i < n / 2; ++i)
        if (!(s[i] == s[n - 1 - i])) return false;
    return true;
}

// ---- 区切り・つなげる ----
inline std::vector<std::string> split(const std::string& s, char sep, bool skip_empty = false) {  // sep で区切る（空の部分も残す。skip_empty で除く）
    std::vector<std::string> r;
    std::size_t b = 0;
    for (std::size_t i = 0; i <= s.size(); ++i)
        if (i == s.size() || s[i] == sep) {
            if (!(skip_empty && i == b)) r.push_back(s.substr(b, i - b));
            b = i + 1;
        }
    return r;
}
template <class V>
std::string join(const std::vector<V>& vs, const std::string& sep = " ") {  // 区切りを挟んでつなげる（数も文字列にする。答えの出力に）
    std::string r;
    for (std::size_t i = 0; i < vs.size(); ++i) {
        if (i) r += sep;
        if constexpr (std::is_same<V, std::string>::value) r += vs[i];
        else if constexpr (std::is_same<V, char>::value) r += vs[i];
        else if constexpr (std::is_integral<V>::value) r += std::to_string(vs[i]);
        else {
            std::ostringstream os;
            os << vs[i];
            r += os.str();
        }
    }
    return r;
}
#endif  // GG_STR_BASE
#ifndef GG_PALINDROME
#define GG_PALINDROME
// ---- 回文（文字列第二世代）  Manacher・区間の回文判定 O(1)・最長の回文・回文の個数・回文にする・回文部分列・回文分割 ----
// 中心の番号 c（0 ≤ c < 2N - 1）：c = 2i は s[i] が中心（奇数長）、c = 2i - 1 は s[i-1] と s[i] の間が中心（偶数長）
namespace str_internal {
inline std::vector<int> manacher_codes(const std::vector<long long>& s) {  // L[c] = 中心 c の最長の回文の長さ
    const int n = (int)s.size();
    if (n == 0) return {};
    std::vector<int> d1(n), d2(n), L(2 * n - 1);
    for (int i = 0, l = 0, r = -1; i < n; ++i) {
        int k = (i > r) ? 1 : std::min(d1[l + r - i], r - i + 1);
        while (i - k >= 0 && i + k < n && s[i - k] == s[i + k]) ++k;
        d1[i] = k--;
        if (i + k > r) l = i - k, r = i + k;
    }
    for (int i = 0, l = 0, r = -1; i < n; ++i) {
        int k = (i > r) ? 0 : std::min(d2[l + r - i + 1], r - i + 1);
        while (i - k - 1 >= 0 && i + k < n && s[i - k - 1] == s[i + k]) ++k;
        d2[i] = k--;
        if (i + k > r) l = i - k - 1, r = i + k;
    }
    for (int i = 0; i < n; ++i) L[2 * i] = 2 * d1[i] - 1;
    for (int i = 1; i < n; ++i) L[2 * i - 1] = 2 * d2[i];
    return L;
}
}  // namespace str_internal

template <class Seq>
std::vector<int> center_palindrome_lengths(const Seq& s) {  // 中心 c（0〜2N-2）ごとの最長の回文の長さ（Manacher）  O(N)
    return str_internal::manacher_codes(str_internal::codes_of(s));
}

struct Palindromes {
    Palindromes() = default;
    template <class Seq>
    explicit Palindromes(const Seq& s) : n_((int)str_internal::seq(s).size()), L_(center_palindrome_lengths(s)) {}
    int size() const { return n_; }
    bool is_palindrome(int l, int r) const {  // s[l, r) が回文か（空・1 文字は回文）  O(1)
        str_internal::check_range("Palindromes::is_palindrome", l, r, n_);
        return r - l <= 1 || L_[l + r - 1] >= r - l;
    }
    int radius_odd(int i) const {  // s[i] が中心の最長の回文の半径 r（長さ 2r - 1）
        if (i < 0 || i >= n_) str_internal::fail("Palindromes::radius_odd", "i は 0〜N-1");
        return (L_[2 * i] + 1) / 2;
    }
    int radius_even(int i) const {  // s[i-1] と s[i] の間が中心の最長の回文の半径 r（長さ 2r。i = 0, N なら 0）
        if (i < 0 || i > n_) str_internal::fail("Palindromes::radius_even", "i は 0〜N");
        return (i == 0 || i == n_) ? 0 : L_[2 * i - 1] / 2;
    }
    int length_at_center(int c) const {  // 中心 c の最長の回文の長さ
        if (c < 0 || c >= 2 * n_ - 1) str_internal::fail("Palindromes::length_at_center", "c は 0〜2N-2");
        return L_[c];
    }
    StrSpan span_at_center(int c) const {  // 中心 c の最長の回文の位置
        const int len = length_at_center(c), l = (c + 1 - len) / 2;
        return {l, l + len};
    }
    StrSpan longest() const {  // 最長の回文の部分文字列（同じ長さなら左）  O(N)
        StrSpan best{0, 0};
        for (int c = 0; c < 2 * n_ - 1; ++c)
            if (L_[c] > best.len() || (L_[c] == best.len() && (c + 1 - L_[c]) / 2 < best.l)) best = span_at_center(c);
        return best;
    }
    long long count_all() const {  // 回文である部分文字列（空でない）の数。位置が違えば別  O(N)
        long long t = 0;
        for (int c = 0; c < 2 * n_ - 1; ++c) t += (L_[c] + (c % 2 == 0)) / 2;
        return t;
    }
    int longest_prefix_palindrome() const {  // s[0, r) が回文となる最大の r  O(N)
        for (int r = n_; r > 0; --r)
            if (is_palindrome(0, r)) return r;
        return 0;
    }
    int longest_suffix_palindrome() const {  // s[N - k, N) が回文となる最大の k  O(N)
        for (int l = 0; l < n_; ++l)
            if (is_palindrome(l, n_)) return n_ - l;
        return 0;
    }
    std::vector<int> longest_starting_at() const {  // r[i] = s[i, i + len) が回文となる最大の len  O(N)
        std::vector<int> best(n_, -1), r(n_);
        for (int c = 0; c < 2 * n_ - 1; ++c) {
            const int l = (c + 1 - L_[c]) / 2;
            best[l] = std::max(best[l], c);
        }
        for (int i = 0, mx = -1; i < n_; ++i) mx = std::max(mx, best[i]), r[i] = mx - 2 * i + 1;
        return r;
    }
    std::vector<int> longest_ending_at() const {  // r[j] = s[j - len + 1, j + 1) が回文となる最大の len（j が最後の文字）  O(N)
        std::vector<int> best(n_, 2 * n_), r(n_);
        for (int c = 0; c < 2 * n_ - 1; ++c) {
            const int e = (c + 1 - L_[c]) / 2 + L_[c] - 1;  // 最後の文字の位置
            best[e] = std::min(best[e], c);
        }
        for (int j = n_ - 1, mn = 2 * n_; j >= 0; --j) mn = std::min(mn, best[j]), r[j] = 2 * j - mn + 1;
        return r;
    }

   private:
    int n_ = 0;
    std::vector<int> L_;
};

// ---- 回文にする・回文の個数 ----
template <class Seq>
str_internal::seq_t<Seq> min_palindrome_with_prefix(const Seq& s0) {  // s で始まる最短の回文（後ろに足す。ABC398 F）  O(N)
    auto s = str_internal::seq_t<Seq>(str_internal::seq(s0));
    const int n = (int)s.size(), k = Palindromes(s).longest_suffix_palindrome();
    for (int i = n - k - 1; i >= 0; --i) s.push_back(s[i]);
    return s;
}
template <class Seq>
str_internal::seq_t<Seq> min_palindrome_with_suffix(const Seq& s0) {  // s で終わる最短の回文（前に足す）  O(N)
    const auto& s = str_internal::seq(s0);
    const int n = (int)s.size(), k = Palindromes(s).longest_prefix_palindrome();
    str_internal::seq_t<Seq> r;
    for (int i = n - 1; i >= k; --i) r.push_back(s[i]);
    r.insert(r.end(), s.begin(), s.end());
    return r;
}
template <class Seq>
int min_changes_to_palindrome(const Seq& s0) {  // 1 文字ずつ書き換えて回文にする最小の回数  O(N)
    const auto& s = str_internal::seq(s0);
    int c = 0;
    for (std::size_t i = 0, n = s.size(); i < n / 2; ++i) c += !(s[i] == s[n - 1 - i]);
    return c;
}
template <class Seq>
bool can_rearrange_to_palindrome(const Seq& s0) {  // 並べ替えて回文にできるか（奇数個の文字が 1 種類以下）  O(N log N)
    auto c = str_internal::codes_of(s0);
    std::sort(c.begin(), c.end());
    int odd = 0;
    for (std::size_t i = 0, j; i < c.size(); i = j) {
        for (j = i; j < c.size() && c[j] == c[i];) ++j;
        odd += (int)((j - i) & 1);
    }
    return odd <= 1;
}
template <class Seq>
int longest_palindromic_subsequence_length(const Seq& s0) {  // 最長回文部分列（飛び飛び）の長さ  O(N^2) 時間・O(N) メモリ
    const auto s = str_internal::codes_of(s0);
    const int n = (int)s.size();
    std::vector<int> nxt(n + 1, 0), cur(n + 1, 0);  // nxt[j] = LPS(s[i+1..j])
    for (int i = n - 1; i >= 0; --i) {
        cur[i] = 1;
        for (int j = i + 1; j < n; ++j) cur[j] = s[i] == s[j] ? (j - 1 >= i + 1 ? nxt[j - 1] : 0) + 2 : std::max(nxt[j], cur[j - 1]);
        std::swap(cur, nxt);
    }
    return n ? nxt[n - 1] : 0;
}
template <class Seq>
str_internal::seq_t<Seq> longest_palindromic_subsequence(const Seq& s0) {  // 最長回文部分列（飛び飛び）を 1 つ  O(N^2) 時間・2N^2 バイト（N ≤ 10^4）
    const auto& s = str_internal::seq(s0);
    const auto c = str_internal::codes_of(s0);
    const int n = (int)c.size();
    if (n > 10000) str_internal::fail("longest_palindromic_subsequence", "N ≤ 10^4 で使う（長さだけなら longest_palindromic_subsequence_length）");
    std::vector<unsigned short> dp((std::size_t)n * n, 0);  // dp[i * n + j] = LPS(s[i..j])（i ≤ j）
    auto at = [&](int i, int j) -> unsigned short& { return dp[(std::size_t)i * n + j]; };
    for (int i = n - 1; i >= 0; --i) {
        at(i, i) = 1;
        for (int j = i + 1; j < n; ++j) at(i, j) = c[i] == c[j] ? (unsigned short)((j - 1 >= i + 1 ? at(i + 1, j - 1) : 0) + 2) : std::max(at(i + 1, j), at(i, j - 1));
    }
    str_internal::seq_t<Seq> left, right;
    int i = 0, j = n - 1;
    while (i <= j) {
        if (i == j) {
            left.push_back(s[i]);
            break;
        }
        if (c[i] == c[j]) left.push_back(s[i]), right.push_back(s[j]), ++i, --j;
        else if (at(i + 1, j) >= at(i, j - 1)) ++i;
        else --j;
    }
    left.insert(left.end(), right.rbegin(), right.rend());
    return left;
}
template <class Seq>
std::vector<StrSpan> palindrome_partition_min(const Seq& s) {  // 回文に分ける最小の個数で分けたもの（.size() が個数）  O(N^2)
    const Palindromes P(s);
    const int n = P.size();
    std::vector<int> dp(n + 1, 1 << 30), from(n + 1, -1);
    dp[0] = 0;
    for (int r = 1; r <= n; ++r)
        for (int l = 0; l < r; ++l)
            if (dp[l] + 1 < dp[r] && P.is_palindrome(l, r)) dp[r] = dp[l] + 1, from[r] = l;
    std::vector<StrSpan> parts;
    for (int r = n; r > 0; r = from[r]) parts.push_back({from[r], r});
    std::reverse(parts.begin(), parts.end());
    return parts;
}
template <class T, class Seq>
T count_palindromic_subsequences(const Seq& s0) {  // 回文である部分列（空でない。位置が違えば別）の数  O(N^2) 時間・O(N) メモリ
    const auto c = str_internal::codes_of(s0);
    const int n = (int)c.size();
    std::vector<T> nxt(n + 1, T(0)), cur(n + 1, T(0));  // f(i, j) = s[i..j] の中の数。nxt は i + 1 の行、添字は j + 1
    for (int i = n - 1; i >= 0; --i) {
        cur[i] = T(0), cur[i + 1] = T(1);  // f(i, i - 1) = 0, f(i, i) = 1
        for (int j = i + 1; j < n; ++j) {
            const T A = nxt[j + 1], B = cur[j], C = nxt[j];  // f(i+1, j), f(i, j-1), f(i+1, j-1)
            cur[j + 1] = c[i] == c[j] ? str_internal::checked_add(str_internal::checked_add(A, B, "count_palindromic_subsequences"), T(1), "count_palindromic_subsequences")
                                      : str_internal::checked_add(A, B - C, "count_palindromic_subsequences");  // B ≥ C なので先に引く（途中で溢れない）
        }
        std::swap(cur, nxt);
    }
    return n ? nxt[n] : T(0);
}
#endif  // GG_PALINDROME

#ifndef GG_STR_BASE
#define GG_STR_BASE
// ---- 文字列基盤（文字列第二世代）  文字の種類 Alphabet・文字の数と位置・ランレングス・変換・比較・区切り ----
// 列（string / vector<整数>）を受け取る関数はどちらでもよい。文字列リテラル（"abc"）も渡せる
// 範囲は半開区間 [l, r)（0 始まり）。位置が「ない」ときは -1。使える範囲の外の引数はメッセージを出して止める
namespace str_internal {
[[noreturn]] inline void fail(const char* func, const std::string& msg) {
    std::cerr << func << ": " << msg << "\n";
    std::abort();
}
// 文字列リテラルを std::string に（それ以外はそのまま参照を返す）
template <class T>
const T& seq(const T& s) { return s; }
template <std::size_t N>
std::string seq(const char (&s)[N]) { return std::string(s); }
inline std::string seq(const char* s) { return std::string(s); }
template <class T>
using seq_t = std::decay_t<decltype(seq(std::declval<const T&>()))>;
template <class T>
using elem_t = std::decay_t<decltype(std::declval<const seq_t<T>&>()[0])>;
// 要素を「値」に（char は unsigned char として 0〜255。負の添字にしない）
template <class T>
long long code_of(T x) {
    if constexpr (std::is_same<T, char>::value || std::is_same<T, signed char>::value || std::is_same<T, unsigned char>::value)
        return (long long)(unsigned char)x;
    else
        return (long long)x;
}
// 範囲の確認（止めるときの文字列作りは別の関数に分けて、確認そのものを軽くする）
[[noreturn]] __attribute__((noinline, cold)) inline void range_fail(const char* func, long long l, long long r, long long n) {
    fail(func, "範囲 [l, r) は 0 ≤ l ≤ r ≤ 長さ で指定する（l = " + std::to_string(l) + ", r = " + std::to_string(r) + ", 長さ = " + std::to_string(n) + "）");
}
inline void check_range(const char* func, long long l, long long r, long long n) {
    if (__builtin_expect(l < 0 || r > n || l > r, 0)) range_fail(func, l, r, n);
}
// a[lo, hi) の中で x 以上の最初の位置（a は昇順）。分岐の少ない二分探索
inline int lower_index(const int* a, int lo, int hi, int x) {
    int len = hi - lo;
    if (len <= 0) return lo;
    const int* b = a + lo;
    while (len > 1) {
        const int half = len >> 1;
        b += (b[half - 1] < x) ? half : 0;
        len -= half;
    }
    return (int)(b - a) + (*b < x);
}
// 昇順の a の中で x より大きい最初の位置（std::upper_bound と同じ。デバッグビルドでも O(log N)）
template <class T>
int upper_index(const std::vector<T>& a, const T& x) {
    int lo = 0, hi = (int)a.size();
    while (lo < hi) {
        const int mid = (int)((unsigned)(lo + hi) >> 1);
        if (!(x < a[mid])) lo = mid + 1;
        else hi = mid;
    }
    return lo;
}
// 数え上げの足し算：T が整数型なら溢れたら止める（mint・double などはそのまま）
template <class T>
T checked_add(const T& a, const T& b, const char* func) {
    if constexpr (std::is_integral<T>::value) {
        T r;
        if (__builtin_add_overflow(a, b, &r)) fail(func, "値が型に収まらない（答えを mod で求めるなら math_v2 の ModInt998 などを使う）");
        return r;
    } else {
        return a + b;
    }
}
// 数え上げの途中の足し算：T が整数型なら溢れたことを印 oa に残すだけ（途中で溢れても答えに使わなければ止めない。最後に印を見て止める）
template <class T>
void add_count(T& a, bool& oa, const T& b, bool ob) {
    if constexpr (std::is_integral<T>::value) {
        T r;
        oa = __builtin_add_overflow(a, b, &r) || oa || ob;
        a = r;
    } else {
        a = a + b;
    }
}
[[noreturn]] inline void count_overflow(const char* func) { fail(func, "値が型に収まらない（答えを mod で求めるなら math_v2 の ModInt998 などを使う）"); }
// 列 → 値の列（char は 0〜255）
template <class Seq>
std::vector<long long> codes_of(const Seq& s0) {
    const auto& s = seq(s0);
    std::vector<long long> c(s.size());
    for (std::size_t i = 0; i < s.size(); ++i) c[i] = code_of(s[i]);
    return c;
}
// 辞書順の比較（-1 / 0 / 1）。char は unsigned char として比べる（std::string と同じ順）
template <class A, class B>
int lex_compare(const A& a, const B& b) {
    const std::size_t n = a.size(), m = b.size(), k = n < m ? n : m;
    for (std::size_t i = 0; i < k; ++i) {
        const long long x = code_of(a[i]), y = code_of(b[i]);
        if (x != y) return x < y ? -1 : 1;
    }
    return n == m ? 0 : (n < m ? -1 : 1);
}
}  // namespace str_internal

// ---- 部分文字列の位置 [l, r) ----
struct StrSpan {
    int l = 0, r = 0;
    int len() const { return r - l; }
    bool empty() const { return r <= l; }
    friend bool operator==(const StrSpan& a, const StrSpan& b) { return a.l == b.l && a.r == b.r; }
    friend bool operator!=(const StrSpan& a, const StrSpan& b) { return !(a == b); }
    friend std::ostream& operator<<(std::ostream& os, const StrSpan& sp) { return os << "[" << sp.l << ", " << sp.r << ")"; }  // デバッグ用
};
template <class Seq>
str_internal::seq_t<Seq> substr_of(const Seq& s0, StrSpan sp) {  // s[sp.l, sp.r) のコピー（string なら string、vector なら vector）
    const auto& s = str_internal::seq(s0);
    str_internal::check_range("substr_of", sp.l, sp.r, (long long)s.size());
    return str_internal::seq_t<Seq>(s.begin() + sp.l, s.begin() + sp.r);
}

// ---- 文字の種類（Alphabet）：文字・値 → 番号 0〜size()-1。番号の順が辞書順として扱われる ----
//   Alphabet::lower()（a〜z。既定）/ upper()（A〜Z）/ digits()（0〜9）/ letters()（A〜Z, a〜z の 52 種。ASCII 順）
//   alnum()（0〜9, A〜Z, a〜z の 62 種。ASCII 順）/ ascii()（0〜127）/ bytes()（0〜255。日本語など 128 以上の char も）
//   of("()[]")（並べた文字だけ。並べた順が番号の順 = 辞書順）/ ints(K)（整数 0〜K-1）/ ints(lo, hi)（整数 [lo, hi)）
//   範囲外の文字・値を id に渡すとメッセージを出して止まる（contains で先に確かめられる）
struct Alphabet {
    Alphabet() : Alphabet(of("abcdefghijklmnopqrstuvwxyz")) {}
    static Alphabet lower() { return of("abcdefghijklmnopqrstuvwxyz"); }
    static Alphabet upper() { return of("ABCDEFGHIJKLMNOPQRSTUVWXYZ"); }
    static Alphabet digits() { return of("0123456789"); }
    static Alphabet letters() { return of("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"); }
    static Alphabet alnum() { return of("0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"); }
    static Alphabet ascii() { return code_range(0, 128); }
    static Alphabet bytes() { return code_range(0, 256); }
    static Alphabet of(const std::string& cs) {  // 並べた文字（同じ文字を 2 回並べると止まる）
        if (cs.empty()) str_internal::fail("Alphabet::of", "文字が 1 つもない");
        Alphabet A(0);
        A.kind_ = CHARS;
        for (std::size_t i = 0; i < cs.size(); ++i) {
            const unsigned char c = (unsigned char)cs[i];
            if (A.table_[c] >= 0) str_internal::fail("Alphabet::of", std::string("同じ文字が 2 回ある: '") + cs[i] + "'");
            A.table_[c] = (short)i;
            A.chars_[i] = c;
        }
        A.size_ = (int)cs.size();
        return A;
    }
    static Alphabet ints(long long K) { return ints(0, K); }  // 整数 0〜K-1
    static Alphabet ints(long long lo, long long hi) {         // 整数 [lo, hi)
        if (hi <= lo || (unsigned long long)hi - (unsigned long long)lo > 1000000000ULL)
            str_internal::fail("Alphabet::ints", "lo < hi かつ hi - lo ≤ 10^9 で指定する（大きな値は compress_values で 0〜K-1 に）");
        Alphabet A(0);
        A.kind_ = INTS;
        A.lo_ = lo;
        A.size_ = (int)(hi - lo);
        return A;
    }

    int size() const { return size_; }
    bool is_chars() const { return kind_ == CHARS; }  // 文字の表で決めたものか（of・lower など）。false なら ints
    template <class T>
    bool contains(T x) const { return index_of(str_internal::code_of(x)) >= 0; }
    template <class T>
    int id(T x) const {  // 文字・値の番号（範囲外は止める）
        const long long v = str_internal::code_of(x);
        const int r = index_of(v);
        if (r < 0) out_of_range(v);
        return r;
    }
    long long value(int i) const {  // 番号 i の値（文字なら文字コード 0〜255）
        check_index(i);
        return kind_ == CHARS ? (long long)chars_[i] : lo_ + i;
    }
    char ch(int i) const {  // 番号 i の文字（ints で 0〜255 の外の値なら止める）
        const long long v = value(i);
        if (v < 0 || v > 255) str_internal::fail("Alphabet::ch", "番号 " + std::to_string(i) + " の値 " + std::to_string(v) + " は文字にできない");
        return (char)(unsigned char)v;
    }
    std::string chars() const {  // 文字を番号の順に並べたもの（ints なら 0〜255 の範囲のときだけ）
        std::string s;
        for (int i = 0; i < size_; ++i) s += ch(i);
        return s;
    }

   private:
    enum Kind : unsigned char { CHARS, INTS };
    explicit Alphabet(int) { table_.fill(-1); }
    static Alphabet code_range(int lo, int hi) {
        Alphabet A(0);
        A.kind_ = CHARS;
        for (int c = lo; c < hi; ++c) A.table_[c] = (short)(c - lo), A.chars_[c - lo] = (unsigned char)c;
        A.size_ = hi - lo;
        return A;
    }
    int index_of(long long v) const {
        if (kind_ == CHARS) return (unsigned long long)v < 256ULL ? table_[(std::size_t)v] : -1;
        const unsigned long long d = (unsigned long long)v - (unsigned long long)lo_;
        return d < (unsigned long long)size_ ? (int)d : -1;
    }
    void check_index(int i) const {
        if (i < 0 || i >= size_) str_internal::fail("Alphabet", "番号 " + std::to_string(i) + " は 0〜" + std::to_string(size_ - 1) + " の外");
    }
    [[noreturn]] void out_of_range(long long v) const {
        std::string what = "値 " + std::to_string(v) + " は";
        if (kind_ == CHARS && v >= 32 && v < 127) what = std::string("文字 '") + (char)v + "'（コード " + std::to_string(v) + "）は";
        str_internal::fail("Alphabet::id", what + "文字の種類に入っていない（種類の数 " + std::to_string(size_) + "）");
    }
    Kind kind_ = CHARS;
    int size_ = 0;
    long long lo_ = 0;
    std::array<short, 256> table_{};
    std::array<unsigned char, 256> chars_{};
};

template <class Seq>
std::vector<int> to_ids(const Seq& s0, const Alphabet& A = Alphabet::lower()) {  // 各文字の番号の列
    const auto& s = str_internal::seq(s0);
    std::vector<int> r(s.size());
    for (std::size_t i = 0; i < s.size(); ++i) r[i] = A.id(s[i]);
    return r;
}
inline std::string from_ids(const std::vector<int>& ids, const Alphabet& A = Alphabet::lower()) {  // 番号の列 → 文字列
    std::string s(ids.size(), '\0');
    for (std::size_t i = 0; i < ids.size(); ++i) s[i] = A.ch(ids[i]);
    return s;
}

// ---- 数の列を 0〜K-1 に圧縮（大小の順を保つ）：Alphabet::ints(K) と組み合わせる ----
struct CompressedValues {
    std::vector<int> ids;            // 各要素の番号（0〜K-1。値の小さい順）
    std::vector<long long> values;   // 番号 → 元の値（昇順。K = values.size()）
    int size() const { return (int)values.size(); }
    Alphabet alphabet() const { return Alphabet::ints(values.empty() ? 1 : (long long)values.size()); }
};
template <class Seq>
CompressedValues compress_values(const Seq& s0) {
    const auto& s = str_internal::seq(s0);
    const std::size_t n = s.size();
    std::vector<std::pair<long long, int>> vi(n);
    for (std::size_t i = 0; i < n; ++i) vi[i] = {str_internal::code_of(s[i]), (int)i};
    std::sort(vi.begin(), vi.end());
    CompressedValues c;
    c.ids.resize(n);
    for (std::size_t i = 0; i < n; ++i) {
        if (i == 0 || vi[i].first != vi[i - 1].first) c.values.push_back(vi[i].first);
        c.ids[vi[i].second] = (int)c.values.size() - 1;
    }
    return c;
}

// ---- 文字の数 ----
template <class Seq>
std::vector<int> count_chars(const Seq& s0, const Alphabet& A = Alphabet::lower()) {  // 番号ごとの個数（長さ A.size()）
    const auto& s = str_internal::seq(s0);
    std::vector<int> c(A.size(), 0);
    for (std::size_t i = 0; i < s.size(); ++i) ++c[A.id(s[i])];
    return c;
}

// ---- 区間の文字の数・次 / 前の出現を O(1) で（前処理 O(Nσ)・メモリ (N+1)σ 個の int。N = 2×10^5, σ = 26 で約 21MB） ----
//   メモリを抑えたいとき・σ が大きいときは CharPositions（O(log N)・メモリ O(N + σ)）
struct CharPrefix {
    CharPrefix() = default;
    template <class Seq>
    explicit CharPrefix(const Seq& s0, const Alphabet& A = Alphabet::lower()) : A_(A) {
        const auto& s = str_internal::seq(s0);
        n_ = (int)s.size(), sg_ = A_.size();
        if ((long long)(n_ + 1) * sg_ > 600000000LL) str_internal::fail("CharPrefix", "(N + 1)σ が 6×10^8 を超える（メモリが足りない。CharPositions を使う）");
        cnt_.assign((std::size_t)(n_ + 1) * sg_, 0);
        std::vector<int> id(n_);
        start_.assign(sg_ + 1, 0);
        for (int i = 0; i < n_; ++i) {
            const int* src = cnt_.data() + (std::size_t)i * sg_;
            int* dst = cnt_.data() + (std::size_t)(i + 1) * sg_;
            std::copy(src, src + sg_, dst);
            id[i] = A_.id(s[i]);
            ++dst[id[i]], ++start_[id[i] + 1];
        }
        for (int k = 0; k < sg_; ++k) start_[k + 1] += start_[k];
        pos_.resize(n_);
        std::vector<int> cur(start_.begin(), start_.end() - 1);
        for (int i = 0; i < n_; ++i) pos_[cur[id[i]]++] = i;
    }
    int size() const { return n_; }
    const Alphabet& alphabet() const { return A_; }
    template <class C>
    int count(int l, int r, C c) const {  // [l, r) にある c の数  O(1)
        str_internal::check_range("CharPrefix::count", l, r, n_);
        const int k = A_.id(c);
        return cnt_[(std::size_t)r * sg_ + k] - cnt_[(std::size_t)l * sg_ + k];
    }
    template <class C>
    bool has(int l, int r, C c) const { return count(l, r, c) > 0; }  // [l, r) に c があるか  O(1)
    std::vector<int> counts(int l, int r) const {  // [l, r) の番号ごとの数  O(σ)
        str_internal::check_range("CharPrefix::counts", l, r, n_);
        std::vector<int> c(sg_);
        for (int k = 0; k < sg_; ++k) c[k] = cnt_[(std::size_t)r * sg_ + k] - cnt_[(std::size_t)l * sg_ + k];
        return c;
    }
    int distinct(int l, int r) const {  // [l, r) にある文字の種類の数  O(σ)
        str_internal::check_range("CharPrefix::distinct", l, r, n_);
        int d = 0;
        for (int k = 0; k < sg_; ++k) d += cnt_[(std::size_t)r * sg_ + k] != cnt_[(std::size_t)l * sg_ + k];
        return d;
    }
    template <class C>
    int next(int i, C c) const {  // i 以上で最初の c の位置（なければ -1。0 ≤ i ≤ N）  O(1)
        check_pos("CharPrefix::next", i);
        const int k = A_.id(c), m = cnt_[(std::size_t)i * sg_ + k];
        return start_[k] + m < start_[k + 1] ? pos_[start_[k] + m] : -1;
    }
    template <class C>
    int prev(int i, C c) const {  // i より前（[0, i)）で最後の c の位置（なければ -1。0 ≤ i ≤ N）  O(1)
        check_pos("CharPrefix::prev", i);
        const int k = A_.id(c), m = cnt_[(std::size_t)i * sg_ + k];
        return m > 0 ? pos_[start_[k] + m - 1] : -1;
    }
    template <class C>
    int kth(C c, long long k) const {  // c の k 番目（0 から）の位置（なければ -1）  O(1)
        const int id = A_.id(c);
        return (k >= 0 && k < start_[id + 1] - start_[id]) ? pos_[start_[id] + k] : -1;
    }
    template <class C>
    int total(C c) const {  // 全体での c の数  O(1)
        const int k = A_.id(c);
        return start_[k + 1] - start_[k];
    }

   private:
    void check_pos(const char* func, int i) const {
        if (i < 0 || i > n_) str_internal::fail(func, "位置 " + std::to_string(i) + " は 0〜N（= " + std::to_string(n_) + "）の外");
    }
    Alphabet A_;
    int n_ = 0, sg_ = 0;
    std::vector<int> cnt_;           // cnt_[i * σ + k] = [0, i) にある番号 k の数
    std::vector<int> start_, pos_;   // 番号ごとの出現位置（pos_[start_[k] .. start_[k+1])）
};

// ---- 文字の出現位置（前処理 O(N + σ)・メモリ O(N + σ)。σ が大きくても使える） ----
struct CharPositions {
    CharPositions() = default;
    template <class Seq>
    explicit CharPositions(const Seq& s0, const Alphabet& A = Alphabet::lower()) : A_(A) {
        const auto& s = str_internal::seq(s0);
        n_ = (int)s.size();
        const int sg = A_.size();
        std::vector<int> id(n_);
        start_.assign(sg + 1, 0);
        for (int i = 0; i < n_; ++i) id[i] = A_.id(s[i]), ++start_[id[i] + 1];
        for (int k = 0; k < sg; ++k) start_[k + 1] += start_[k];
        pos_.resize(n_), occ_.resize(n_);
        std::vector<int> cur(start_.begin(), start_.end() - 1);
        for (int i = 0; i < n_; ++i) occ_[i] = cur[id[i]] - start_[id[i]], pos_[cur[id[i]]++] = i;
    }
    int size() const { return n_; }
    const Alphabet& alphabet() const { return A_; }
    template <class C>
    int next(int i, C c) const {  // i 以上で最初の c の位置（なければ -1。0 ≤ i ≤ N）  O(log N)
        check_pos("CharPositions::next", i);
        const int k = A_.id(c), j = str_internal::lower_index(pos_.data(), start_[k], start_[k + 1], i);
        return j < start_[k + 1] ? pos_[j] : -1;
    }
    template <class C>
    int prev(int i, C c) const {  // i より前（[0, i)）で最後の c の位置（なければ -1。0 ≤ i ≤ N）  O(log N)
        check_pos("CharPositions::prev", i);
        const int k = A_.id(c), j = str_internal::lower_index(pos_.data(), start_[k], start_[k + 1], i);
        return j > start_[k] ? pos_[j - 1] : -1;
    }
    template <class C>
    int kth(C c, long long k) const {  // c の k 番目（0 から）の位置（なければ -1）  O(1)
        const int id = A_.id(c);
        return (k >= 0 && k < start_[id + 1] - start_[id]) ? pos_[start_[id] + k] : -1;
    }
    template <class C>
    int count(int l, int r, C c) const {  // [l, r) にある c の数  O(log N)
        str_internal::check_range("CharPositions::count", l, r, n_);
        const int k = A_.id(c);
        return str_internal::lower_index(pos_.data(), start_[k], start_[k + 1], r) - str_internal::lower_index(pos_.data(), start_[k], start_[k + 1], l);
    }
    template <class C>
    bool has(int l, int r, C c) const { return count(l, r, c) > 0; }
    template <class C>
    int total(C c) const {  // 全体での c の数  O(1)
        const int k = A_.id(c);
        return start_[k + 1] - start_[k];
    }
    template <class C>
    std::vector<int> positions(C c) const {  // c の位置の一覧（昇順）
        const int k = A_.id(c);
        return std::vector<int>(pos_.begin() + start_[k], pos_.begin() + start_[k + 1]);
    }
    int occurrence_index(int i) const {  // s[i] が「その文字の何番目（0 から）の出現」か  O(1)
        if (i < 0 || i >= n_) str_internal::fail("CharPositions::occurrence_index", "位置 " + std::to_string(i) + " は 0〜N-1 の外");
        return occ_[i];
    }

   private:
    void check_pos(const char* func, int i) const {
        if (i < 0 || i > n_) str_internal::fail(func, "位置 " + std::to_string(i) + " は 0〜N（= " + std::to_string(n_) + "）の外");
    }
    Alphabet A_;
    int n_ = 0;
    std::vector<int> start_, pos_, occ_;
};

// ---- ランレングス圧縮：連続する同じ値を {value, len, start} にまとめる ----
template <class T>
struct RunBlock {
    T value{};
    long long len = 0;    // 続く長さ
    long long start = 0;  // 元の列での開始位置
    friend std::ostream& operator<<(std::ostream& os, const RunBlock& b) { return os << "{" << b.value << ", " << b.len << ", " << b.start << "}"; }  // デバッグ用
};
template <class Seq>
std::vector<RunBlock<str_internal::elem_t<Seq>>> run_length(const Seq& s0) {
    const auto& s = str_internal::seq(s0);
    std::vector<RunBlock<str_internal::elem_t<Seq>>> r;
    std::size_t cnt = s.empty() ? 0 : 1;
    for (std::size_t i = 1; i < s.size(); ++i) cnt += !(s[i] == s[i - 1]);
    r.reserve(cnt);
    for (std::size_t i = 0, n = s.size(); i < n;) {
        std::size_t j = i + 1;
        while (j < n && s[j] == s[i]) ++j;
        r.push_back({s[i], (long long)(j - i), (long long)i});
        i = j;
    }
    return r;
}
template <class T>
auto run_length_decode(const std::vector<RunBlock<T>>& runs) {  // 元に戻す（T = char なら string、それ以外は vector<T>）
    using R = std::conditional_t<std::is_same<T, char>::value, std::string, std::vector<T>>;
    long long total = 0;
    for (const auto& b : runs) {
        if (b.len < 0) str_internal::fail("run_length_decode", "len が負");
        total += b.len;
        if (total > 2000000000LL) str_internal::fail("run_length_decode", "元の長さが 2×10^9 を超える");
    }
    R r;
    r.reserve((std::size_t)total);
    for (const auto& b : runs) r.insert(r.end(), (std::size_t)b.len, b.value);
    return r;
}

// ---- 変換 ----
inline std::string to_lower(std::string s) {  // 英大文字を小文字に（ほかはそのまま）
    for (char& c : s)
        if ('A' <= c && c <= 'Z') c = (char)(c - 'A' + 'a');
    return s;
}
inline std::string to_upper(std::string s) {
    for (char& c : s)
        if ('a' <= c && c <= 'z') c = (char)(c - 'a' + 'A');
    return s;
}
inline std::string swap_case(std::string s) {  // 大文字と小文字を入れ替える
    for (char& c : s) {
        if ('A' <= c && c <= 'Z') c = (char)(c - 'A' + 'a');
        else if ('a' <= c && c <= 'z') c = (char)(c - 'a' + 'A');
    }
    return s;
}
template <class Seq>
str_internal::seq_t<Seq> caesar_shift(const Seq& s0, long long k, const Alphabet& A = Alphabet::lower()) {  // A の文字を番号で k ずらす（負も。A にない文字はそのまま）
    auto s = str_internal::seq_t<Seq>(str_internal::seq(s0));
    const long long sg = A.size(), sh = ((k % sg) + sg) % sg;
    for (auto& c : s)
        if (A.contains(c)) c = (std::decay_t<decltype(c)>)A.value((int)((A.id(c) + sh) % sg));
    return s;
}
// 文字の置き換えを重ねる：replace(c, d) で「今 c になっている文字を全部 d に」。apply(s) で s に当てはめる（ABC342 C）
struct CharMap {
    explicit CharMap(const Alphabet& A = Alphabet::lower()) : A_(A), f_(A.size()) {
        for (int i = 0; i < A_.size(); ++i) f_[i] = i;
    }
    template <class C, class D>
    void replace(C c, D d) {  // 今 c になっているものを全部 d に  O(σ)
        const int x = A_.id(c), y = A_.id(d);
        if (x == y) return;
        for (int& v : f_)
            if (v == x) v = y;
    }
    template <class C, class D>
    void set(C c, D d) { f_[A_.id(c)] = A_.id(d); }  // 元の文字 c の行き先を d にする（ほかは変えない）
    template <class C>
    auto operator()(C c) const {  // 元の文字 c の行き先（A にない文字はそのまま）
        return A_.contains(c) ? (C)A_.value(f_[A_.id(c)]) : c;
    }
    template <class Seq>
    str_internal::seq_t<Seq> apply(const Seq& s0) const {  // s の各文字を行き先に（A にない文字はそのまま）
        auto s = str_internal::seq_t<Seq>(str_internal::seq(s0));
        for (auto& c : s) c = (*this)(c);
        return s;
    }
    CharMap compose(const CharMap& g) const {  // 「この置き換えのあとに g」を 1 つにしたもの
        CharMap h(A_);
        for (int i = 0; i < A_.size(); ++i) h.f_[i] = g.f_[f_[i]];
        return h;
    }
    const std::vector<int>& table() const { return f_; }  // 番号 → 行き先の番号

   private:
    Alphabet A_;
    std::vector<int> f_;
};
template <class Seq>
str_internal::seq_t<Seq> reversed(const Seq& s0) {  // 反転したコピー
    const auto& s = str_internal::seq(s0);
    return str_internal::seq_t<Seq>(s.rbegin(), s.rend());
}
template <class Seq>
str_internal::seq_t<Seq> min_with_reverse(const Seq& s0) {  // s と反転のうち辞書順で小さい方（向きを区別しないときの代表。ABC310 C）
    auto s = str_internal::seq_t<Seq>(str_internal::seq(s0));
    auto t = reversed(s);
    return str_internal::lex_compare(t, s) < 0 ? t : s;
}
template <class Seq>
str_internal::seq_t<Seq> sorted_chars(const Seq& s0) {  // 並べ替えたもの（並べ替えて同じかの比較用）
    auto s = str_internal::seq_t<Seq>(str_internal::seq(s0));
    std::sort(s.begin(), s.end(), [](const auto& a, const auto& b) { return str_internal::code_of(a) < str_internal::code_of(b); });
    return s;
}
template <class Seq>
str_internal::seq_t<Seq> repeat_string(const Seq& s0, long long k) {  // s を k 回つなげる
    const auto& s = str_internal::seq(s0);
    if (k < 0) str_internal::fail("repeat_string", "k ≥ 0 で使う");
    if (!s.empty() && k > 2000000000LL / (long long)s.size()) str_internal::fail("repeat_string", "長さが 2×10^9 を超える");
    str_internal::seq_t<Seq> r;
    r.reserve((std::size_t)(k * (long long)s.size()));
    for (long long i = 0; i < k; ++i) r.insert(r.end(), s.begin(), s.end());
    return r;
}

// ---- 比較 ----
template <class A, class B>
int common_prefix_length(const A& a0, const B& b0) {  // 共通の接頭辞の長さ
    const auto& a = str_internal::seq(a0);
    const auto& b = str_internal::seq(b0);
    const std::size_t k = a.size() < b.size() ? a.size() : b.size();
    std::size_t i = 0;
    while (i < k && str_internal::code_of(a[i]) == str_internal::code_of(b[i])) ++i;
    return (int)i;
}
template <class A, class B>
int common_suffix_length(const A& a0, const B& b0) {  // 共通の接尾辞の長さ
    const auto& a = str_internal::seq(a0);
    const auto& b = str_internal::seq(b0);
    const std::size_t k = a.size() < b.size() ? a.size() : b.size();
    std::size_t i = 0;
    while (i < k && str_internal::code_of(a[a.size() - 1 - i]) == str_internal::code_of(b[b.size() - 1 - i])) ++i;
    return (int)i;
}
template <class A, class B>
int hamming_distance(const A& a0, const B& b0) {  // 同じ位置で違う文字の数（長さが同じこと）
    const auto& a = str_internal::seq(a0);
    const auto& b = str_internal::seq(b0);
    if (a.size() != b.size()) str_internal::fail("hamming_distance", "長さが違う");
    int d = 0;
    for (std::size_t i = 0; i < a.size(); ++i) d += str_internal::code_of(a[i]) != str_internal::code_of(b[i]);
    return d;
}
template <class T, class S>
bool is_subsequence(const T& t0, const S& s0) {  // t が s の部分列（飛び飛びでよい）か  O(|s|)。何度も聞くなら 06_subsequence の NextTable
    const auto& t = str_internal::seq(t0);
    const auto& s = str_internal::seq(s0);
    std::size_t j = 0;
    for (std::size_t i = 0; i < s.size() && j < t.size(); ++i)
        if (str_internal::code_of(s[i]) == str_internal::code_of(t[j])) ++j;
    return j == t.size();
}
template <class Seq>
bool is_palindrome(const Seq& s0) {  // 全体が回文か  O(N)。区間を何度も聞くなら 04_palindrome の Palindromes
    const auto& s = str_internal::seq(s0);
    for (std::size_t i = 0, n = s.size(); i < n / 2; ++i)
        if (!(s[i] == s[n - 1 - i])) return false;
    return true;
}

// ---- 区切り・つなげる ----
inline std::vector<std::string> split(const std::string& s, char sep, bool skip_empty = false) {  // sep で区切る（空の部分も残す。skip_empty で除く）
    std::vector<std::string> r;
    std::size_t b = 0;
    for (std::size_t i = 0; i <= s.size(); ++i)
        if (i == s.size() || s[i] == sep) {
            if (!(skip_empty && i == b)) r.push_back(s.substr(b, i - b));
            b = i + 1;
        }
    return r;
}
template <class V>
std::string join(const std::vector<V>& vs, const std::string& sep = " ") {  // 区切りを挟んでつなげる（数も文字列にする。答えの出力に）
    std::string r;
    for (std::size_t i = 0; i < vs.size(); ++i) {
        if (i) r += sep;
        if constexpr (std::is_same<V, std::string>::value) r += vs[i];
        else if constexpr (std::is_same<V, char>::value) r += vs[i];
        else if constexpr (std::is_integral<V>::value) r += std::to_string(vs[i]);
        else {
            std::ostringstream os;
            os << vs[i];
            r += os.str();
        }
    }
    return r;
}
#endif  // GG_STR_BASE
#ifndef GG_LEX_ORDER
#define GG_LEX_ORDER
// ---- 辞書順（文字列第二世代）  つなげて最小・独自の文字の順・最小の回転・列名（a, …, z, aa）・辞書の番号・次の文字列 ----
// 辞書順は「先頭から比べ、最初に違う文字で決まる。片方がもう片方の接頭辞なら短い方が小さい」（std::string と同じ）
namespace str_internal {
__extension__ typedef __int128 lex_i128;
// a + b と b + a の辞書順（-1 / 0 / 1）を、つなげた列を作らずに
template <class S>
int concat_cmp(const S& a, const S& b) {
    const std::size_t n = a.size(), m = b.size(), t = n + m;
    for (std::size_t i = 0; i < t; ++i) {
        const long long x = code_of(i < n ? a[i] : b[i - n]), y = code_of(i < m ? b[i] : a[i - m]);
        if (x != y) return x < y ? -1 : 1;
    }
    return 0;
}
inline lex_i128 lex_big() { return (lex_i128)1000000000000000000LL * 10; }  // 10^19（long long の最大より大きい）
inline lex_i128 lex_cap(lex_i128 x) { return x > lex_big() ? lex_big() : x; }
// 長さ 0〜t の文字列の数 1 + σ + … + σ^t（10^19 で止める）
inline lex_i128 lex_count_upto(long long sigma, long long t) {
    if (t < 0) return 0;
    if (sigma == 1) return lex_cap((lex_i128)t + 1);
    lex_i128 sum = 0, p = 1;
    for (long long u = 0; u <= t; ++u) {
        sum = lex_cap(sum + p);
        p = lex_cap(p * sigma);
        if (sum >= lex_big()) break;
    }
    return sum;
}
}  // namespace str_internal

// ---- つなげて辞書順最小 ----
template <class S>
std::vector<int> concat_order(const std::vector<S>& vs) {  // 全部をつなげて辞書順最小になる並べ方（添字。a + b < b + a の順。同じなら元の順。空の列は最後）  O(N log N × 最長の長さ)
    std::vector<int> idx, empty;
    for (int i = 0; i < (int)vs.size(); ++i) (str_internal::seq(vs[i]).size() ? idx : empty).push_back(i);  // 空の列は順序を壊すので外す
    std::stable_sort(idx.begin(), idx.end(), [&](int i, int j) { return str_internal::concat_cmp(vs[i], vs[j]) < 0; });
    idx.insert(idx.end(), empty.begin(), empty.end());
    return idx;
}
template <class S>
S min_concat(const std::vector<S>& vs) {  // 全部をつなげた辞書順最小の列
    S r;
    for (int i : concat_order(vs)) r.insert(r.end(), vs[i].begin(), vs[i].end());
    return r;
}

// ---- 独自の文字の順での辞書順（order に並べた順。order にない文字は止める） ----
struct LexLess {
    explicit LexLess(const Alphabet& A) : A_(A) {}
    template <class X, class Y>
    int compare(const X& a0, const Y& b0) const {  // -1 / 0 / 1
        const auto& a = str_internal::seq(a0);
        const auto& b = str_internal::seq(b0);
        const std::size_t k = std::min(a.size(), b.size());
        for (std::size_t i = 0; i < k; ++i) {
            const int x = A_.id(a[i]), y = A_.id(b[i]);
            if (x != y) return x < y ? -1 : 1;
        }
        return a.size() == b.size() ? 0 : (a.size() < b.size() ? -1 : 1);
    }
    template <class X, class Y>
    bool operator()(const X& a, const Y& b) const { return compare(a, b) < 0; }  // sort の比較関数に

   private:
    Alphabet A_;
};
inline LexLess lex_less_by(const std::string& order) { return LexLess(Alphabet::of(order)); }  // sort(v.begin(), v.end(), lex_less_by("bac…"))
inline LexLess lex_less_by(const Alphabet& A) { return LexLess(A); }

// ---- 回転 ----
namespace str_internal {
template <class Seq>
int extreme_rotation(const Seq& s0, bool want_min) {
    const auto s = codes_of(s0);
    const int n = (int)s.size();
    if (n == 0) return 0;
    int i = 0, j = 1, k = 0;
    while (i < n && j < n && k < n) {
        const int ia = i + k < n ? i + k : i + k - n, ib = j + k < n ? j + k : j + k - n;
        const long long a = s[ia], b = s[ib];
        if (a == b) {
            ++k;
            continue;
        }
        if ((a > b) == want_min) i += k + 1;
        else j += k + 1;
        if (i == j) ++j;
        k = 0;
    }
    int r = std::min(i, j);
    // 周期があれば最小の位置に（s = t^m なら r mod |t|）
    std::vector<int> pi(n);
    for (int x = 1; x < n; ++x) {
        int q = pi[x - 1];
        while (q > 0 && s[x] != s[q]) q = pi[q - 1];
        pi[x] = q + (s[x] == s[q]);
    }
    const int p = n - pi[n - 1];
    if (n % p == 0) r %= p;
    return r;
}
}  // namespace str_internal
template <class Seq>
int min_rotation(const Seq& s) { return str_internal::extreme_rotation(s, true); }  // 辞書順最小の回転 s[k..] + s[..k] の k（同じなら最小の k）  O(N)
template <class Seq>
int max_rotation(const Seq& s) { return str_internal::extreme_rotation(s, false); }  // 辞書順最大の回転の k  O(N)

// ---- 列名（a, b, …, z, aa, ab, …）：全単射の記数法 ----
inline std::string bijective_to_string(long long n, const Alphabet& A = Alphabet::lower()) {  // 1 → "a"、26 → "z"、27 → "aa"（0 なら空）。ABC171 C
    if (n < 0) str_internal::fail("bijective_to_string", "n ≥ 0 で使う");
    const long long k = A.size();
    std::string r;
    while (n > 0) {
        --n;
        r += A.ch((int)(n % k));
        n /= k;
    }
    std::reverse(r.begin(), r.end());
    return r;
}
template <class Seq>
long long bijective_from_string(const Seq& s0, const Alphabet& A = Alphabet::lower()) {  // "a" → 1、"z" → 26、"aa" → 27（ABC285 C）。long long に収まらなければ止める
    const auto& s = str_internal::seq(s0);
    const long long k = A.size();
    str_internal::lex_i128 v = 0;
    for (std::size_t i = 0; i < s.size(); ++i) {
        v = v * k + (A.id(s[i]) + 1);
        if (v > (str_internal::lex_i128)std::numeric_limits<long long>::max()) str_internal::fail("bijective_from_string", "long long に収まらない");
    }
    return (long long)v;
}

// ---- 辞書の番号：長さ 1〜max_len のすべての文字列（A の文字）を辞書順に並べたときの番号（0 から） ----
template <class Seq>
long long dict_rank(const Seq& s0, long long max_len, const Alphabet& A = Alphabet::lower()) {  // s より小さい文字列の数  O(|s| + max_len)
    const auto& s = str_internal::seq(s0);
    const long long m = (long long)s.size(), sg = A.size();
    if (m < 1 || m > max_len) str_internal::fail("dict_rank", "1 ≤ |s| ≤ max_len で使う");
    str_internal::lex_i128 r = m - 1;  // s の真の接頭辞（空でない）
    for (long long d = 0; d < m; ++d) {
        const int c = A.id(s[d]);
        if (c) r += (str_internal::lex_i128)c * str_internal::lex_count_upto(sg, max_len - d - 1);
        if (r > (str_internal::lex_i128)std::numeric_limits<long long>::max()) str_internal::fail("dict_rank", "番号が long long に収まらない");
    }
    return (long long)r;
}
inline std::string dict_unrank(long long k, long long max_len, const Alphabet& A = Alphabet::lower()) {  // 番号 k（0 から）の文字列  O(max_len)
    if (k < 0 || max_len < 1) str_internal::fail("dict_unrank", "k ≥ 0, max_len ≥ 1 で使う");
    const long long sg = A.size();
    if (str_internal::lex_count_upto(sg, max_len) - 1 <= k) str_internal::fail("dict_unrank", "k が文字列の数以上");
    std::string r;
    str_internal::lex_i128 x = k;
    for (long long d = 0; d < max_len; ++d) {
        const str_internal::lex_i128 block = str_internal::lex_count_upto(sg, max_len - d - 1);  // r + c で始まる文字列の数
        const long long c = (long long)(x / block);
        r += A.ch((int)c);
        x -= (str_internal::lex_i128)c * block;
        if (x == 0) break;
        --x;  // r + c そのものを飛ばす
    }
    return r;
}
template <class Seq>
str_internal::seq_t<Seq> next_string(const Seq& s0, const Alphabet& A = Alphabet::lower()) {  // 同じ長さで辞書順の次（最後なら空）  O(N)
    auto s = str_internal::seq_t<Seq>(str_internal::seq(s0));
    for (int i = (int)s.size() - 1; i >= 0; --i) {
        const int c = A.id(s[i]);
        if (c + 1 < A.size()) {
            s[i] = (typename str_internal::seq_t<Seq>::value_type)A.value(c + 1);
            for (int j = i + 1; j < (int)s.size(); ++j) s[j] = (typename str_internal::seq_t<Seq>::value_type)A.value(0);
            return s;
        }
    }
    return {};
}
#endif  // GG_LEX_ORDER

#ifndef GG_STR_BASE
#define GG_STR_BASE
// ---- 文字列基盤（文字列第二世代）  文字の種類 Alphabet・文字の数と位置・ランレングス・変換・比較・区切り ----
// 列（string / vector<整数>）を受け取る関数はどちらでもよい。文字列リテラル（"abc"）も渡せる
// 範囲は半開区間 [l, r)（0 始まり）。位置が「ない」ときは -1。使える範囲の外の引数はメッセージを出して止める
namespace str_internal {
[[noreturn]] inline void fail(const char* func, const std::string& msg) {
    std::cerr << func << ": " << msg << "\n";
    std::abort();
}
// 文字列リテラルを std::string に（それ以外はそのまま参照を返す）
template <class T>
const T& seq(const T& s) { return s; }
template <std::size_t N>
std::string seq(const char (&s)[N]) { return std::string(s); }
inline std::string seq(const char* s) { return std::string(s); }
template <class T>
using seq_t = std::decay_t<decltype(seq(std::declval<const T&>()))>;
template <class T>
using elem_t = std::decay_t<decltype(std::declval<const seq_t<T>&>()[0])>;
// 要素を「値」に（char は unsigned char として 0〜255。負の添字にしない）
template <class T>
long long code_of(T x) {
    if constexpr (std::is_same<T, char>::value || std::is_same<T, signed char>::value || std::is_same<T, unsigned char>::value)
        return (long long)(unsigned char)x;
    else
        return (long long)x;
}
// 範囲の確認（止めるときの文字列作りは別の関数に分けて、確認そのものを軽くする）
[[noreturn]] __attribute__((noinline, cold)) inline void range_fail(const char* func, long long l, long long r, long long n) {
    fail(func, "範囲 [l, r) は 0 ≤ l ≤ r ≤ 長さ で指定する（l = " + std::to_string(l) + ", r = " + std::to_string(r) + ", 長さ = " + std::to_string(n) + "）");
}
inline void check_range(const char* func, long long l, long long r, long long n) {
    if (__builtin_expect(l < 0 || r > n || l > r, 0)) range_fail(func, l, r, n);
}
// a[lo, hi) の中で x 以上の最初の位置（a は昇順）。分岐の少ない二分探索
inline int lower_index(const int* a, int lo, int hi, int x) {
    int len = hi - lo;
    if (len <= 0) return lo;
    const int* b = a + lo;
    while (len > 1) {
        const int half = len >> 1;
        b += (b[half - 1] < x) ? half : 0;
        len -= half;
    }
    return (int)(b - a) + (*b < x);
}
// 昇順の a の中で x より大きい最初の位置（std::upper_bound と同じ。デバッグビルドでも O(log N)）
template <class T>
int upper_index(const std::vector<T>& a, const T& x) {
    int lo = 0, hi = (int)a.size();
    while (lo < hi) {
        const int mid = (int)((unsigned)(lo + hi) >> 1);
        if (!(x < a[mid])) lo = mid + 1;
        else hi = mid;
    }
    return lo;
}
// 数え上げの足し算：T が整数型なら溢れたら止める（mint・double などはそのまま）
template <class T>
T checked_add(const T& a, const T& b, const char* func) {
    if constexpr (std::is_integral<T>::value) {
        T r;
        if (__builtin_add_overflow(a, b, &r)) fail(func, "値が型に収まらない（答えを mod で求めるなら math_v2 の ModInt998 などを使う）");
        return r;
    } else {
        return a + b;
    }
}
// 数え上げの途中の足し算：T が整数型なら溢れたことを印 oa に残すだけ（途中で溢れても答えに使わなければ止めない。最後に印を見て止める）
template <class T>
void add_count(T& a, bool& oa, const T& b, bool ob) {
    if constexpr (std::is_integral<T>::value) {
        T r;
        oa = __builtin_add_overflow(a, b, &r) || oa || ob;
        a = r;
    } else {
        a = a + b;
    }
}
[[noreturn]] inline void count_overflow(const char* func) { fail(func, "値が型に収まらない（答えを mod で求めるなら math_v2 の ModInt998 などを使う）"); }
// 列 → 値の列（char は 0〜255）
template <class Seq>
std::vector<long long> codes_of(const Seq& s0) {
    const auto& s = seq(s0);
    std::vector<long long> c(s.size());
    for (std::size_t i = 0; i < s.size(); ++i) c[i] = code_of(s[i]);
    return c;
}
// 辞書順の比較（-1 / 0 / 1）。char は unsigned char として比べる（std::string と同じ順）
template <class A, class B>
int lex_compare(const A& a, const B& b) {
    const std::size_t n = a.size(), m = b.size(), k = n < m ? n : m;
    for (std::size_t i = 0; i < k; ++i) {
        const long long x = code_of(a[i]), y = code_of(b[i]);
        if (x != y) return x < y ? -1 : 1;
    }
    return n == m ? 0 : (n < m ? -1 : 1);
}
}  // namespace str_internal

// ---- 部分文字列の位置 [l, r) ----
struct StrSpan {
    int l = 0, r = 0;
    int len() const { return r - l; }
    bool empty() const { return r <= l; }
    friend bool operator==(const StrSpan& a, const StrSpan& b) { return a.l == b.l && a.r == b.r; }
    friend bool operator!=(const StrSpan& a, const StrSpan& b) { return !(a == b); }
    friend std::ostream& operator<<(std::ostream& os, const StrSpan& sp) { return os << "[" << sp.l << ", " << sp.r << ")"; }  // デバッグ用
};
template <class Seq>
str_internal::seq_t<Seq> substr_of(const Seq& s0, StrSpan sp) {  // s[sp.l, sp.r) のコピー（string なら string、vector なら vector）
    const auto& s = str_internal::seq(s0);
    str_internal::check_range("substr_of", sp.l, sp.r, (long long)s.size());
    return str_internal::seq_t<Seq>(s.begin() + sp.l, s.begin() + sp.r);
}

// ---- 文字の種類（Alphabet）：文字・値 → 番号 0〜size()-1。番号の順が辞書順として扱われる ----
//   Alphabet::lower()（a〜z。既定）/ upper()（A〜Z）/ digits()（0〜9）/ letters()（A〜Z, a〜z の 52 種。ASCII 順）
//   alnum()（0〜9, A〜Z, a〜z の 62 種。ASCII 順）/ ascii()（0〜127）/ bytes()（0〜255。日本語など 128 以上の char も）
//   of("()[]")（並べた文字だけ。並べた順が番号の順 = 辞書順）/ ints(K)（整数 0〜K-1）/ ints(lo, hi)（整数 [lo, hi)）
//   範囲外の文字・値を id に渡すとメッセージを出して止まる（contains で先に確かめられる）
struct Alphabet {
    Alphabet() : Alphabet(of("abcdefghijklmnopqrstuvwxyz")) {}
    static Alphabet lower() { return of("abcdefghijklmnopqrstuvwxyz"); }
    static Alphabet upper() { return of("ABCDEFGHIJKLMNOPQRSTUVWXYZ"); }
    static Alphabet digits() { return of("0123456789"); }
    static Alphabet letters() { return of("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"); }
    static Alphabet alnum() { return of("0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"); }
    static Alphabet ascii() { return code_range(0, 128); }
    static Alphabet bytes() { return code_range(0, 256); }
    static Alphabet of(const std::string& cs) {  // 並べた文字（同じ文字を 2 回並べると止まる）
        if (cs.empty()) str_internal::fail("Alphabet::of", "文字が 1 つもない");
        Alphabet A(0);
        A.kind_ = CHARS;
        for (std::size_t i = 0; i < cs.size(); ++i) {
            const unsigned char c = (unsigned char)cs[i];
            if (A.table_[c] >= 0) str_internal::fail("Alphabet::of", std::string("同じ文字が 2 回ある: '") + cs[i] + "'");
            A.table_[c] = (short)i;
            A.chars_[i] = c;
        }
        A.size_ = (int)cs.size();
        return A;
    }
    static Alphabet ints(long long K) { return ints(0, K); }  // 整数 0〜K-1
    static Alphabet ints(long long lo, long long hi) {         // 整数 [lo, hi)
        if (hi <= lo || (unsigned long long)hi - (unsigned long long)lo > 1000000000ULL)
            str_internal::fail("Alphabet::ints", "lo < hi かつ hi - lo ≤ 10^9 で指定する（大きな値は compress_values で 0〜K-1 に）");
        Alphabet A(0);
        A.kind_ = INTS;
        A.lo_ = lo;
        A.size_ = (int)(hi - lo);
        return A;
    }

    int size() const { return size_; }
    bool is_chars() const { return kind_ == CHARS; }  // 文字の表で決めたものか（of・lower など）。false なら ints
    template <class T>
    bool contains(T x) const { return index_of(str_internal::code_of(x)) >= 0; }
    template <class T>
    int id(T x) const {  // 文字・値の番号（範囲外は止める）
        const long long v = str_internal::code_of(x);
        const int r = index_of(v);
        if (r < 0) out_of_range(v);
        return r;
    }
    long long value(int i) const {  // 番号 i の値（文字なら文字コード 0〜255）
        check_index(i);
        return kind_ == CHARS ? (long long)chars_[i] : lo_ + i;
    }
    char ch(int i) const {  // 番号 i の文字（ints で 0〜255 の外の値なら止める）
        const long long v = value(i);
        if (v < 0 || v > 255) str_internal::fail("Alphabet::ch", "番号 " + std::to_string(i) + " の値 " + std::to_string(v) + " は文字にできない");
        return (char)(unsigned char)v;
    }
    std::string chars() const {  // 文字を番号の順に並べたもの（ints なら 0〜255 の範囲のときだけ）
        std::string s;
        for (int i = 0; i < size_; ++i) s += ch(i);
        return s;
    }

   private:
    enum Kind : unsigned char { CHARS, INTS };
    explicit Alphabet(int) { table_.fill(-1); }
    static Alphabet code_range(int lo, int hi) {
        Alphabet A(0);
        A.kind_ = CHARS;
        for (int c = lo; c < hi; ++c) A.table_[c] = (short)(c - lo), A.chars_[c - lo] = (unsigned char)c;
        A.size_ = hi - lo;
        return A;
    }
    int index_of(long long v) const {
        if (kind_ == CHARS) return (unsigned long long)v < 256ULL ? table_[(std::size_t)v] : -1;
        const unsigned long long d = (unsigned long long)v - (unsigned long long)lo_;
        return d < (unsigned long long)size_ ? (int)d : -1;
    }
    void check_index(int i) const {
        if (i < 0 || i >= size_) str_internal::fail("Alphabet", "番号 " + std::to_string(i) + " は 0〜" + std::to_string(size_ - 1) + " の外");
    }
    [[noreturn]] void out_of_range(long long v) const {
        std::string what = "値 " + std::to_string(v) + " は";
        if (kind_ == CHARS && v >= 32 && v < 127) what = std::string("文字 '") + (char)v + "'（コード " + std::to_string(v) + "）は";
        str_internal::fail("Alphabet::id", what + "文字の種類に入っていない（種類の数 " + std::to_string(size_) + "）");
    }
    Kind kind_ = CHARS;
    int size_ = 0;
    long long lo_ = 0;
    std::array<short, 256> table_{};
    std::array<unsigned char, 256> chars_{};
};

template <class Seq>
std::vector<int> to_ids(const Seq& s0, const Alphabet& A = Alphabet::lower()) {  // 各文字の番号の列
    const auto& s = str_internal::seq(s0);
    std::vector<int> r(s.size());
    for (std::size_t i = 0; i < s.size(); ++i) r[i] = A.id(s[i]);
    return r;
}
inline std::string from_ids(const std::vector<int>& ids, const Alphabet& A = Alphabet::lower()) {  // 番号の列 → 文字列
    std::string s(ids.size(), '\0');
    for (std::size_t i = 0; i < ids.size(); ++i) s[i] = A.ch(ids[i]);
    return s;
}

// ---- 数の列を 0〜K-1 に圧縮（大小の順を保つ）：Alphabet::ints(K) と組み合わせる ----
struct CompressedValues {
    std::vector<int> ids;            // 各要素の番号（0〜K-1。値の小さい順）
    std::vector<long long> values;   // 番号 → 元の値（昇順。K = values.size()）
    int size() const { return (int)values.size(); }
    Alphabet alphabet() const { return Alphabet::ints(values.empty() ? 1 : (long long)values.size()); }
};
template <class Seq>
CompressedValues compress_values(const Seq& s0) {
    const auto& s = str_internal::seq(s0);
    const std::size_t n = s.size();
    std::vector<std::pair<long long, int>> vi(n);
    for (std::size_t i = 0; i < n; ++i) vi[i] = {str_internal::code_of(s[i]), (int)i};
    std::sort(vi.begin(), vi.end());
    CompressedValues c;
    c.ids.resize(n);
    for (std::size_t i = 0; i < n; ++i) {
        if (i == 0 || vi[i].first != vi[i - 1].first) c.values.push_back(vi[i].first);
        c.ids[vi[i].second] = (int)c.values.size() - 1;
    }
    return c;
}

// ---- 文字の数 ----
template <class Seq>
std::vector<int> count_chars(const Seq& s0, const Alphabet& A = Alphabet::lower()) {  // 番号ごとの個数（長さ A.size()）
    const auto& s = str_internal::seq(s0);
    std::vector<int> c(A.size(), 0);
    for (std::size_t i = 0; i < s.size(); ++i) ++c[A.id(s[i])];
    return c;
}

// ---- 区間の文字の数・次 / 前の出現を O(1) で（前処理 O(Nσ)・メモリ (N+1)σ 個の int。N = 2×10^5, σ = 26 で約 21MB） ----
//   メモリを抑えたいとき・σ が大きいときは CharPositions（O(log N)・メモリ O(N + σ)）
struct CharPrefix {
    CharPrefix() = default;
    template <class Seq>
    explicit CharPrefix(const Seq& s0, const Alphabet& A = Alphabet::lower()) : A_(A) {
        const auto& s = str_internal::seq(s0);
        n_ = (int)s.size(), sg_ = A_.size();
        if ((long long)(n_ + 1) * sg_ > 600000000LL) str_internal::fail("CharPrefix", "(N + 1)σ が 6×10^8 を超える（メモリが足りない。CharPositions を使う）");
        cnt_.assign((std::size_t)(n_ + 1) * sg_, 0);
        std::vector<int> id(n_);
        start_.assign(sg_ + 1, 0);
        for (int i = 0; i < n_; ++i) {
            const int* src = cnt_.data() + (std::size_t)i * sg_;
            int* dst = cnt_.data() + (std::size_t)(i + 1) * sg_;
            std::copy(src, src + sg_, dst);
            id[i] = A_.id(s[i]);
            ++dst[id[i]], ++start_[id[i] + 1];
        }
        for (int k = 0; k < sg_; ++k) start_[k + 1] += start_[k];
        pos_.resize(n_);
        std::vector<int> cur(start_.begin(), start_.end() - 1);
        for (int i = 0; i < n_; ++i) pos_[cur[id[i]]++] = i;
    }
    int size() const { return n_; }
    const Alphabet& alphabet() const { return A_; }
    template <class C>
    int count(int l, int r, C c) const {  // [l, r) にある c の数  O(1)
        str_internal::check_range("CharPrefix::count", l, r, n_);
        const int k = A_.id(c);
        return cnt_[(std::size_t)r * sg_ + k] - cnt_[(std::size_t)l * sg_ + k];
    }
    template <class C>
    bool has(int l, int r, C c) const { return count(l, r, c) > 0; }  // [l, r) に c があるか  O(1)
    std::vector<int> counts(int l, int r) const {  // [l, r) の番号ごとの数  O(σ)
        str_internal::check_range("CharPrefix::counts", l, r, n_);
        std::vector<int> c(sg_);
        for (int k = 0; k < sg_; ++k) c[k] = cnt_[(std::size_t)r * sg_ + k] - cnt_[(std::size_t)l * sg_ + k];
        return c;
    }
    int distinct(int l, int r) const {  // [l, r) にある文字の種類の数  O(σ)
        str_internal::check_range("CharPrefix::distinct", l, r, n_);
        int d = 0;
        for (int k = 0; k < sg_; ++k) d += cnt_[(std::size_t)r * sg_ + k] != cnt_[(std::size_t)l * sg_ + k];
        return d;
    }
    template <class C>
    int next(int i, C c) const {  // i 以上で最初の c の位置（なければ -1。0 ≤ i ≤ N）  O(1)
        check_pos("CharPrefix::next", i);
        const int k = A_.id(c), m = cnt_[(std::size_t)i * sg_ + k];
        return start_[k] + m < start_[k + 1] ? pos_[start_[k] + m] : -1;
    }
    template <class C>
    int prev(int i, C c) const {  // i より前（[0, i)）で最後の c の位置（なければ -1。0 ≤ i ≤ N）  O(1)
        check_pos("CharPrefix::prev", i);
        const int k = A_.id(c), m = cnt_[(std::size_t)i * sg_ + k];
        return m > 0 ? pos_[start_[k] + m - 1] : -1;
    }
    template <class C>
    int kth(C c, long long k) const {  // c の k 番目（0 から）の位置（なければ -1）  O(1)
        const int id = A_.id(c);
        return (k >= 0 && k < start_[id + 1] - start_[id]) ? pos_[start_[id] + k] : -1;
    }
    template <class C>
    int total(C c) const {  // 全体での c の数  O(1)
        const int k = A_.id(c);
        return start_[k + 1] - start_[k];
    }

   private:
    void check_pos(const char* func, int i) const {
        if (i < 0 || i > n_) str_internal::fail(func, "位置 " + std::to_string(i) + " は 0〜N（= " + std::to_string(n_) + "）の外");
    }
    Alphabet A_;
    int n_ = 0, sg_ = 0;
    std::vector<int> cnt_;           // cnt_[i * σ + k] = [0, i) にある番号 k の数
    std::vector<int> start_, pos_;   // 番号ごとの出現位置（pos_[start_[k] .. start_[k+1])）
};

// ---- 文字の出現位置（前処理 O(N + σ)・メモリ O(N + σ)。σ が大きくても使える） ----
struct CharPositions {
    CharPositions() = default;
    template <class Seq>
    explicit CharPositions(const Seq& s0, const Alphabet& A = Alphabet::lower()) : A_(A) {
        const auto& s = str_internal::seq(s0);
        n_ = (int)s.size();
        const int sg = A_.size();
        std::vector<int> id(n_);
        start_.assign(sg + 1, 0);
        for (int i = 0; i < n_; ++i) id[i] = A_.id(s[i]), ++start_[id[i] + 1];
        for (int k = 0; k < sg; ++k) start_[k + 1] += start_[k];
        pos_.resize(n_), occ_.resize(n_);
        std::vector<int> cur(start_.begin(), start_.end() - 1);
        for (int i = 0; i < n_; ++i) occ_[i] = cur[id[i]] - start_[id[i]], pos_[cur[id[i]]++] = i;
    }
    int size() const { return n_; }
    const Alphabet& alphabet() const { return A_; }
    template <class C>
    int next(int i, C c) const {  // i 以上で最初の c の位置（なければ -1。0 ≤ i ≤ N）  O(log N)
        check_pos("CharPositions::next", i);
        const int k = A_.id(c), j = str_internal::lower_index(pos_.data(), start_[k], start_[k + 1], i);
        return j < start_[k + 1] ? pos_[j] : -1;
    }
    template <class C>
    int prev(int i, C c) const {  // i より前（[0, i)）で最後の c の位置（なければ -1。0 ≤ i ≤ N）  O(log N)
        check_pos("CharPositions::prev", i);
        const int k = A_.id(c), j = str_internal::lower_index(pos_.data(), start_[k], start_[k + 1], i);
        return j > start_[k] ? pos_[j - 1] : -1;
    }
    template <class C>
    int kth(C c, long long k) const {  // c の k 番目（0 から）の位置（なければ -1）  O(1)
        const int id = A_.id(c);
        return (k >= 0 && k < start_[id + 1] - start_[id]) ? pos_[start_[id] + k] : -1;
    }
    template <class C>
    int count(int l, int r, C c) const {  // [l, r) にある c の数  O(log N)
        str_internal::check_range("CharPositions::count", l, r, n_);
        const int k = A_.id(c);
        return str_internal::lower_index(pos_.data(), start_[k], start_[k + 1], r) - str_internal::lower_index(pos_.data(), start_[k], start_[k + 1], l);
    }
    template <class C>
    bool has(int l, int r, C c) const { return count(l, r, c) > 0; }
    template <class C>
    int total(C c) const {  // 全体での c の数  O(1)
        const int k = A_.id(c);
        return start_[k + 1] - start_[k];
    }
    template <class C>
    std::vector<int> positions(C c) const {  // c の位置の一覧（昇順）
        const int k = A_.id(c);
        return std::vector<int>(pos_.begin() + start_[k], pos_.begin() + start_[k + 1]);
    }
    int occurrence_index(int i) const {  // s[i] が「その文字の何番目（0 から）の出現」か  O(1)
        if (i < 0 || i >= n_) str_internal::fail("CharPositions::occurrence_index", "位置 " + std::to_string(i) + " は 0〜N-1 の外");
        return occ_[i];
    }

   private:
    void check_pos(const char* func, int i) const {
        if (i < 0 || i > n_) str_internal::fail(func, "位置 " + std::to_string(i) + " は 0〜N（= " + std::to_string(n_) + "）の外");
    }
    Alphabet A_;
    int n_ = 0;
    std::vector<int> start_, pos_, occ_;
};

// ---- ランレングス圧縮：連続する同じ値を {value, len, start} にまとめる ----
template <class T>
struct RunBlock {
    T value{};
    long long len = 0;    // 続く長さ
    long long start = 0;  // 元の列での開始位置
    friend std::ostream& operator<<(std::ostream& os, const RunBlock& b) { return os << "{" << b.value << ", " << b.len << ", " << b.start << "}"; }  // デバッグ用
};
template <class Seq>
std::vector<RunBlock<str_internal::elem_t<Seq>>> run_length(const Seq& s0) {
    const auto& s = str_internal::seq(s0);
    std::vector<RunBlock<str_internal::elem_t<Seq>>> r;
    std::size_t cnt = s.empty() ? 0 : 1;
    for (std::size_t i = 1; i < s.size(); ++i) cnt += !(s[i] == s[i - 1]);
    r.reserve(cnt);
    for (std::size_t i = 0, n = s.size(); i < n;) {
        std::size_t j = i + 1;
        while (j < n && s[j] == s[i]) ++j;
        r.push_back({s[i], (long long)(j - i), (long long)i});
        i = j;
    }
    return r;
}
template <class T>
auto run_length_decode(const std::vector<RunBlock<T>>& runs) {  // 元に戻す（T = char なら string、それ以外は vector<T>）
    using R = std::conditional_t<std::is_same<T, char>::value, std::string, std::vector<T>>;
    long long total = 0;
    for (const auto& b : runs) {
        if (b.len < 0) str_internal::fail("run_length_decode", "len が負");
        total += b.len;
        if (total > 2000000000LL) str_internal::fail("run_length_decode", "元の長さが 2×10^9 を超える");
    }
    R r;
    r.reserve((std::size_t)total);
    for (const auto& b : runs) r.insert(r.end(), (std::size_t)b.len, b.value);
    return r;
}

// ---- 変換 ----
inline std::string to_lower(std::string s) {  // 英大文字を小文字に（ほかはそのまま）
    for (char& c : s)
        if ('A' <= c && c <= 'Z') c = (char)(c - 'A' + 'a');
    return s;
}
inline std::string to_upper(std::string s) {
    for (char& c : s)
        if ('a' <= c && c <= 'z') c = (char)(c - 'a' + 'A');
    return s;
}
inline std::string swap_case(std::string s) {  // 大文字と小文字を入れ替える
    for (char& c : s) {
        if ('A' <= c && c <= 'Z') c = (char)(c - 'A' + 'a');
        else if ('a' <= c && c <= 'z') c = (char)(c - 'a' + 'A');
    }
    return s;
}
template <class Seq>
str_internal::seq_t<Seq> caesar_shift(const Seq& s0, long long k, const Alphabet& A = Alphabet::lower()) {  // A の文字を番号で k ずらす（負も。A にない文字はそのまま）
    auto s = str_internal::seq_t<Seq>(str_internal::seq(s0));
    const long long sg = A.size(), sh = ((k % sg) + sg) % sg;
    for (auto& c : s)
        if (A.contains(c)) c = (std::decay_t<decltype(c)>)A.value((int)((A.id(c) + sh) % sg));
    return s;
}
// 文字の置き換えを重ねる：replace(c, d) で「今 c になっている文字を全部 d に」。apply(s) で s に当てはめる（ABC342 C）
struct CharMap {
    explicit CharMap(const Alphabet& A = Alphabet::lower()) : A_(A), f_(A.size()) {
        for (int i = 0; i < A_.size(); ++i) f_[i] = i;
    }
    template <class C, class D>
    void replace(C c, D d) {  // 今 c になっているものを全部 d に  O(σ)
        const int x = A_.id(c), y = A_.id(d);
        if (x == y) return;
        for (int& v : f_)
            if (v == x) v = y;
    }
    template <class C, class D>
    void set(C c, D d) { f_[A_.id(c)] = A_.id(d); }  // 元の文字 c の行き先を d にする（ほかは変えない）
    template <class C>
    auto operator()(C c) const {  // 元の文字 c の行き先（A にない文字はそのまま）
        return A_.contains(c) ? (C)A_.value(f_[A_.id(c)]) : c;
    }
    template <class Seq>
    str_internal::seq_t<Seq> apply(const Seq& s0) const {  // s の各文字を行き先に（A にない文字はそのまま）
        auto s = str_internal::seq_t<Seq>(str_internal::seq(s0));
        for (auto& c : s) c = (*this)(c);
        return s;
    }
    CharMap compose(const CharMap& g) const {  // 「この置き換えのあとに g」を 1 つにしたもの
        CharMap h(A_);
        for (int i = 0; i < A_.size(); ++i) h.f_[i] = g.f_[f_[i]];
        return h;
    }
    const std::vector<int>& table() const { return f_; }  // 番号 → 行き先の番号

   private:
    Alphabet A_;
    std::vector<int> f_;
};
template <class Seq>
str_internal::seq_t<Seq> reversed(const Seq& s0) {  // 反転したコピー
    const auto& s = str_internal::seq(s0);
    return str_internal::seq_t<Seq>(s.rbegin(), s.rend());
}
template <class Seq>
str_internal::seq_t<Seq> min_with_reverse(const Seq& s0) {  // s と反転のうち辞書順で小さい方（向きを区別しないときの代表。ABC310 C）
    auto s = str_internal::seq_t<Seq>(str_internal::seq(s0));
    auto t = reversed(s);
    return str_internal::lex_compare(t, s) < 0 ? t : s;
}
template <class Seq>
str_internal::seq_t<Seq> sorted_chars(const Seq& s0) {  // 並べ替えたもの（並べ替えて同じかの比較用）
    auto s = str_internal::seq_t<Seq>(str_internal::seq(s0));
    std::sort(s.begin(), s.end(), [](const auto& a, const auto& b) { return str_internal::code_of(a) < str_internal::code_of(b); });
    return s;
}
template <class Seq>
str_internal::seq_t<Seq> repeat_string(const Seq& s0, long long k) {  // s を k 回つなげる
    const auto& s = str_internal::seq(s0);
    if (k < 0) str_internal::fail("repeat_string", "k ≥ 0 で使う");
    if (!s.empty() && k > 2000000000LL / (long long)s.size()) str_internal::fail("repeat_string", "長さが 2×10^9 を超える");
    str_internal::seq_t<Seq> r;
    r.reserve((std::size_t)(k * (long long)s.size()));
    for (long long i = 0; i < k; ++i) r.insert(r.end(), s.begin(), s.end());
    return r;
}

// ---- 比較 ----
template <class A, class B>
int common_prefix_length(const A& a0, const B& b0) {  // 共通の接頭辞の長さ
    const auto& a = str_internal::seq(a0);
    const auto& b = str_internal::seq(b0);
    const std::size_t k = a.size() < b.size() ? a.size() : b.size();
    std::size_t i = 0;
    while (i < k && str_internal::code_of(a[i]) == str_internal::code_of(b[i])) ++i;
    return (int)i;
}
template <class A, class B>
int common_suffix_length(const A& a0, const B& b0) {  // 共通の接尾辞の長さ
    const auto& a = str_internal::seq(a0);
    const auto& b = str_internal::seq(b0);
    const std::size_t k = a.size() < b.size() ? a.size() : b.size();
    std::size_t i = 0;
    while (i < k && str_internal::code_of(a[a.size() - 1 - i]) == str_internal::code_of(b[b.size() - 1 - i])) ++i;
    return (int)i;
}
template <class A, class B>
int hamming_distance(const A& a0, const B& b0) {  // 同じ位置で違う文字の数（長さが同じこと）
    const auto& a = str_internal::seq(a0);
    const auto& b = str_internal::seq(b0);
    if (a.size() != b.size()) str_internal::fail("hamming_distance", "長さが違う");
    int d = 0;
    for (std::size_t i = 0; i < a.size(); ++i) d += str_internal::code_of(a[i]) != str_internal::code_of(b[i]);
    return d;
}
template <class T, class S>
bool is_subsequence(const T& t0, const S& s0) {  // t が s の部分列（飛び飛びでよい）か  O(|s|)。何度も聞くなら 06_subsequence の NextTable
    const auto& t = str_internal::seq(t0);
    const auto& s = str_internal::seq(s0);
    std::size_t j = 0;
    for (std::size_t i = 0; i < s.size() && j < t.size(); ++i)
        if (str_internal::code_of(s[i]) == str_internal::code_of(t[j])) ++j;
    return j == t.size();
}
template <class Seq>
bool is_palindrome(const Seq& s0) {  // 全体が回文か  O(N)。区間を何度も聞くなら 04_palindrome の Palindromes
    const auto& s = str_internal::seq(s0);
    for (std::size_t i = 0, n = s.size(); i < n / 2; ++i)
        if (!(s[i] == s[n - 1 - i])) return false;
    return true;
}

// ---- 区切り・つなげる ----
inline std::vector<std::string> split(const std::string& s, char sep, bool skip_empty = false) {  // sep で区切る（空の部分も残す。skip_empty で除く）
    std::vector<std::string> r;
    std::size_t b = 0;
    for (std::size_t i = 0; i <= s.size(); ++i)
        if (i == s.size() || s[i] == sep) {
            if (!(skip_empty && i == b)) r.push_back(s.substr(b, i - b));
            b = i + 1;
        }
    return r;
}
template <class V>
std::string join(const std::vector<V>& vs, const std::string& sep = " ") {  // 区切りを挟んでつなげる（数も文字列にする。答えの出力に）
    std::string r;
    for (std::size_t i = 0; i < vs.size(); ++i) {
        if (i) r += sep;
        if constexpr (std::is_same<V, std::string>::value) r += vs[i];
        else if constexpr (std::is_same<V, char>::value) r += vs[i];
        else if constexpr (std::is_integral<V>::value) r += std::to_string(vs[i]);
        else {
            std::ostringstream os;
            os << vs[i];
            r += os.str();
        }
    }
    return r;
}
#endif  // GG_STR_BASE
#ifndef GG_LCS_EDIT
#define GG_LCS_EDIT
// ---- 最長共通部分列（LCS）と編集距離（文字列第二世代） ----
// LCS は飛び飛びに選んだ「部分列」。連続する「部分文字列」の最長共通は 03_suffix_array の longest_common_substring
namespace str_internal {
// LCS のビット並列（Hyyrö）：b の位置をビットにして a を 1 文字ずつ読む。rows に各行の V を残す（復元用）
struct LcsBits {
    int n = 0, m = 0, W = 0;
    std::vector<unsigned long long> V;  // 最後の行（keep_rows なら (n+1) 行ぶん）
    template <class A, class B>
    LcsBits(const A& a0, const B& b0, bool keep_rows, const char* func = "lcs_restore") {
        const auto a = codes_of(a0), b = codes_of(b0);
        n = (int)a.size(), m = (int)b.size(), W = (m + 63) / 64;
        if (keep_rows && (long long)(n + 1) * W > 30000000LL) fail(func, "|a| × |b| が大きすぎる（約 1.9×10^9 まで）。長さだけなら lcs_length");
        // b の値ごとの位置のビット
        std::vector<long long> vals(b);
        std::sort(vals.begin(), vals.end());
        vals.erase(std::unique(vals.begin(), vals.end()), vals.end());
        std::vector<unsigned long long> M(vals.size() * (std::size_t)W, 0);
        auto vid = [&](long long x) {
            const int k = upper_index(vals, x) - 1;
            return (k >= 0 && vals[k] == x) ? k : -1;
        };
        for (int j = 0; j < m; ++j) M[(std::size_t)vid(b[j]) * W + j / 64] |= 1ULL << (j % 64);
        const std::size_t rows = keep_rows ? (std::size_t)n + 1 : 1;
        V.assign(rows * W, ~0ULL);
        std::vector<unsigned long long> cur(W, ~0ULL);
        for (int i = 0; i < n; ++i) {
            const int k = vid(a[i]);
            if (k >= 0) {
                const unsigned long long* mk = M.data() + (std::size_t)k * W;
                unsigned long long carry = 0;
                for (int w = 0; w < W; ++w) {
                    // V' = (V + U) | (V - U)、U = V & M（U は V の一部なので V - U = V & ~M。足し算だけ繰り上がりを運ぶ）
                    const unsigned long long v = cur[w], u = v & mk[w];
                    const unsigned long long s1 = v + u, c1 = s1 < v;
                    const unsigned long long s2 = s1 + carry, c2 = s2 < s1;
                    carry = c1 | c2;
                    cur[w] = s2 | (v & ~mk[w]);
                }
            }
            if (keep_rows) std::copy(cur.begin(), cur.end(), V.begin() + (std::size_t)(i + 1) * W);
        }
        if (!keep_rows) std::copy(cur.begin(), cur.end(), V.begin());
    }
    // LCS(a[0, i), b[0, j)) = j - （行 i の下 j ビットの 1 の数）
    int value(int i, int j) const {
        const unsigned long long* r = V.data() + (std::size_t)i * W;
        int ones = 0;
        for (int w = 0; w < j / 64; ++w) ones += __builtin_popcountll(r[w]);
        if (j % 64) ones += __builtin_popcountll(r[j / 64] & ((1ULL << (j % 64)) - 1));
        return j - ones;
    }
};
// LCS の対応する位置の組（a の位置, b の位置）を前から
template <class A, class B>
std::vector<std::pair<int, int>> lcs_pairs(const A& a0, const B& b0, const char* func) {
    const auto a = codes_of(a0), b = codes_of(b0);
    const LcsBits L(a, b, true, func);
    std::vector<std::pair<int, int>> r;
    int i = L.n, j = L.m, cur = L.value(i, j);
    while (i > 0 && j > 0 && cur > 0) {
        if (a[i - 1] == b[j - 1]) {
            r.push_back({i - 1, j - 1}), --i, --j, --cur;
        } else if (L.value(i - 1, j) == cur) {
            --i;
        } else {
            --j;
        }
    }
    std::reverse(r.begin(), r.end());
    return r;
}
}  // namespace str_internal

template <class A, class B>
int lcs_length(const A& a, const B& b) {  // 最長共通部分列（飛び飛び）の長さ。ビット並列  O(|a||b|/64)
    const str_internal::LcsBits L(a, b, false);
    return L.value(0, L.m);
}
template <class A, class B>
str_internal::seq_t<A> lcs_restore(const A& a0, const B& b) {  // 最長共通部分列を 1 つ（a の型で）。|a| × |b| が約 1.9×10^9 まで（両方 4×10^4 くらい）  O(|a||b|/64)
    const auto& a = str_internal::seq(a0);
    str_internal::seq_t<A> r;
    for (auto [i, j] : str_internal::lcs_pairs(a, b, "lcs_restore")) r.push_back(a[i]), (void)j;
    return r;
}
template <class A, class B>
str_internal::seq_t<A> shortest_common_supersequence(const A& a0, const B& b0) {  // a と b の両方を部分列に含む最短の列（長さ |a| + |b| - LCS）
    const auto& a = str_internal::seq(a0);
    const auto& b = str_internal::seq(b0);
    str_internal::seq_t<A> r;
    int pa = 0, pb = 0;
    using E = typename str_internal::seq_t<A>::value_type;
    auto conv = [&](long long v) {  // b の値を a の要素の型に（収まらなければ止める）
        const E e = (E)v;
        if (str_internal::code_of(e) != v) str_internal::fail("shortest_common_supersequence", "b の値 " + std::to_string(v) + " が a の要素の型に収まらない（a を数の列にする）");
        return e;
    };
    for (auto [i, j] : str_internal::lcs_pairs(a, b, "shortest_common_supersequence")) {
        while (pa < i) r.push_back(a[pa++]);
        while (pb < j) r.push_back(conv(str_internal::code_of(b[pb++])));
        r.push_back(a[i]), ++pa, ++pb;
    }
    while (pa < (int)a.size()) r.push_back(a[pa++]);
    while (pb < (int)b.size()) r.push_back(conv(str_internal::code_of(b[pb++])));
    return r;
}

// ---- 編集距離（a を b にする。挿入 ins・削除 del・置換 sub の費用（0 以上）。既定はすべて 1） ----
namespace str_internal {
constexpr long long ed_inf = std::numeric_limits<long long>::max();
inline long long ed_add(long long x, long long y) { return x > ed_inf - y ? ed_inf : x + y; }  // 0 以上どうしの和（上で止める）
inline void ed_check(const char* func, long long ins, long long del, long long sub) {
    if (ins < 0 || del < 0 || sub < 0) fail(func, "費用 ins・del・sub は 0 以上で使う");
}
}  // namespace str_internal
template <class A, class B>
long long edit_distance(const A& a0, const B& b0, long long ins = 1, long long del = 1, long long sub = 1) {  // O(|a||b|) 時間・O(|b|) メモリ
    using str_internal::ed_add;
    str_internal::ed_check("edit_distance", ins, del, sub);
    const auto a = str_internal::codes_of(a0), b = str_internal::codes_of(b0);
    const int n = (int)a.size(), m = (int)b.size();
    std::vector<long long> dp(m + 1);
    dp[0] = 0;
    for (int j = 1; j <= m; ++j) dp[j] = ed_add(dp[j - 1], ins);
    for (int i = 1; i <= n; ++i) {
        long long diag = dp[0];
        dp[0] = ed_add(dp[0], del);
        for (int j = 1; j <= m; ++j) {
            const long long up = dp[j];
            dp[j] = std::min({ed_add(up, del), ed_add(dp[j - 1], ins), a[i - 1] == b[j - 1] ? diag : ed_add(diag, sub)});
            diag = up;
        }
    }
    if (dp[m] == str_internal::ed_inf) str_internal::fail("edit_distance", "答えが long long に収まらない");
    return dp[m];
}
struct EditOp {
    char type;  // '=' そのまま（a[i] == b[j]）/ 'S' a[i] を b[j] に置換 / 'D' a[i] を削除 / 'I' b[j] を挿入（a[i] の前に）
    int i, j;   // a・b での位置（D の j・I の i は「今いる位置」）
};
template <class A, class B>
std::vector<EditOp> edit_script(const A& a0, const B& b0, long long ins = 1, long long del = 1, long long sub = 1) {  // 最小の費用の操作の列（前から）  O(|a||b|) 時間・メモリ
    using str_internal::ed_add;
    str_internal::ed_check("edit_script", ins, del, sub);
    const auto a = str_internal::codes_of(a0), b = str_internal::codes_of(b0);
    const int n = (int)a.size(), m = (int)b.size();
    if ((long long)(n + 1) * (m + 1) > 200000000LL) str_internal::fail("edit_script", "|a| × |b| が 2×10^8 を超える（距離だけなら edit_distance）");
    std::vector<unsigned char> from((std::size_t)(n + 1) * (m + 1));  // 0 = 対角, 1 = 上（削除）, 2 = 左（挿入）
    std::vector<long long> prev(m + 1), cur(m + 1);
    prev[0] = 0, from[0] = 2;
    for (int j = 1; j <= m; ++j) prev[j] = ed_add(prev[j - 1], ins), from[j] = 2;
    for (int i = 1; i <= n; ++i) {
        cur[0] = ed_add(prev[0], del), from[(std::size_t)i * (m + 1)] = 1;
        for (int j = 1; j <= m; ++j) {
            long long best = a[i - 1] == b[j - 1] ? prev[j - 1] : ed_add(prev[j - 1], sub);
            unsigned char f = 0;
            if (ed_add(prev[j], del) < best) best = ed_add(prev[j], del), f = 1;
            if (ed_add(cur[j - 1], ins) < best) best = ed_add(cur[j - 1], ins), f = 2;
            cur[j] = best, from[(std::size_t)i * (m + 1) + j] = f;
        }
        std::swap(prev, cur);
    }
    if (prev[m] == str_internal::ed_inf) str_internal::fail("edit_script", "費用の合計が long long に収まらない");
    std::vector<EditOp> ops;
    for (int i = n, j = m; i > 0 || j > 0;) {
        const unsigned char f = from[(std::size_t)i * (m + 1) + j];
        if (i > 0 && j > 0 && f == 0) ops.push_back({a[i - 1] == b[j - 1] ? '=' : 'S', i - 1, j - 1}), --i, --j;
        else if (i > 0 && (j == 0 || f == 1)) ops.push_back({'D', i - 1, j}), --i;
        else ops.push_back({'I', i, j - 1}), --j;
    }
    std::reverse(ops.begin(), ops.end());
    return ops;
}
template <class A, class B>
int edit_distance_bounded(const A& a0, const B& b0, int k) {  // すべて費用 1 の編集距離。k を超えるなら -1  O((|a| + |b|) k)
    const auto a = str_internal::codes_of(a0), b = str_internal::codes_of(b0);
    const int n = (int)a.size(), m = (int)b.size();
    if (k < 0) str_internal::fail("edit_distance_bounded", "k ≥ 0 で使う");
    if (std::abs(n - m) > k) return -1;
    k = std::min(k, std::max(n, m));  // 距離は max(|a|, |b|) 以下なので、それより大きい k は同じ
    const int W = 2 * k + 1, over = k + 1;  // over = 「k を超えた」
    std::vector<int> prev(W, over), cur(W, over);  // 列 j = i + (d - k)
    for (int d = 0; d < W; ++d) {
        const int j = d - k;
        if (j >= 0 && j <= m) prev[d] = j;
    }
    for (int i = 1; i <= n; ++i) {
        for (int d = 0; d < W; ++d) {
            const int j = i + d - k;
            int v = over;
            if (j >= 0 && j <= m) {
                if (j == 0) v = i;
                else {
                    v = prev[d] + (a[i - 1] == b[j - 1] ? 0 : 1);         // (i-1, j-1) は前の行の同じ d
                    if (d + 1 < W) v = std::min(v, prev[d + 1] + 1);      // (i-1, j)：削除
                    if (d > 0) v = std::min(v, cur[d - 1] + 1);           // (i, j-1)：挿入
                }
            }
            cur[d] = std::min(v, over);
        }
        std::swap(prev, cur);
    }
    const int r = prev[m - n + k];
    return r <= k ? r : -1;
}
template <class A, class B>
bool within_one_edit(const A& a0, const B& b0) {  // 挿入・削除・置換を 1 回以下で一致するか  O(N)
    const auto a = str_internal::codes_of(a0), b = str_internal::codes_of(b0);
    const int n = (int)a.size(), m = (int)b.size();
    if (std::abs(n - m) > 1) return false;
    int i = 0;
    while (i < n && i < m && a[i] == b[i]) ++i;
    if (n == m) {
        int j = n - 1;
        while (j > i && a[j] == b[j]) --j;
        return j <= i;
    }
    const auto& L = n > m ? a : b;
    const auto& S = n > m ? b : a;
    for (int j = i; j < (int)S.size(); ++j)
        if (L[j + 1] != S[j]) return false;
    return true;
}
#endif  // GG_LCS_EDIT

#ifndef GG_STR_BASE
#define GG_STR_BASE
// ---- 文字列基盤（文字列第二世代）  文字の種類 Alphabet・文字の数と位置・ランレングス・変換・比較・区切り ----
// 列（string / vector<整数>）を受け取る関数はどちらでもよい。文字列リテラル（"abc"）も渡せる
// 範囲は半開区間 [l, r)（0 始まり）。位置が「ない」ときは -1。使える範囲の外の引数はメッセージを出して止める
namespace str_internal {
[[noreturn]] inline void fail(const char* func, const std::string& msg) {
    std::cerr << func << ": " << msg << "\n";
    std::abort();
}
// 文字列リテラルを std::string に（それ以外はそのまま参照を返す）
template <class T>
const T& seq(const T& s) { return s; }
template <std::size_t N>
std::string seq(const char (&s)[N]) { return std::string(s); }
inline std::string seq(const char* s) { return std::string(s); }
template <class T>
using seq_t = std::decay_t<decltype(seq(std::declval<const T&>()))>;
template <class T>
using elem_t = std::decay_t<decltype(std::declval<const seq_t<T>&>()[0])>;
// 要素を「値」に（char は unsigned char として 0〜255。負の添字にしない）
template <class T>
long long code_of(T x) {
    if constexpr (std::is_same<T, char>::value || std::is_same<T, signed char>::value || std::is_same<T, unsigned char>::value)
        return (long long)(unsigned char)x;
    else
        return (long long)x;
}
// 範囲の確認（止めるときの文字列作りは別の関数に分けて、確認そのものを軽くする）
[[noreturn]] __attribute__((noinline, cold)) inline void range_fail(const char* func, long long l, long long r, long long n) {
    fail(func, "範囲 [l, r) は 0 ≤ l ≤ r ≤ 長さ で指定する（l = " + std::to_string(l) + ", r = " + std::to_string(r) + ", 長さ = " + std::to_string(n) + "）");
}
inline void check_range(const char* func, long long l, long long r, long long n) {
    if (__builtin_expect(l < 0 || r > n || l > r, 0)) range_fail(func, l, r, n);
}
// a[lo, hi) の中で x 以上の最初の位置（a は昇順）。分岐の少ない二分探索
inline int lower_index(const int* a, int lo, int hi, int x) {
    int len = hi - lo;
    if (len <= 0) return lo;
    const int* b = a + lo;
    while (len > 1) {
        const int half = len >> 1;
        b += (b[half - 1] < x) ? half : 0;
        len -= half;
    }
    return (int)(b - a) + (*b < x);
}
// 昇順の a の中で x より大きい最初の位置（std::upper_bound と同じ。デバッグビルドでも O(log N)）
template <class T>
int upper_index(const std::vector<T>& a, const T& x) {
    int lo = 0, hi = (int)a.size();
    while (lo < hi) {
        const int mid = (int)((unsigned)(lo + hi) >> 1);
        if (!(x < a[mid])) lo = mid + 1;
        else hi = mid;
    }
    return lo;
}
// 数え上げの足し算：T が整数型なら溢れたら止める（mint・double などはそのまま）
template <class T>
T checked_add(const T& a, const T& b, const char* func) {
    if constexpr (std::is_integral<T>::value) {
        T r;
        if (__builtin_add_overflow(a, b, &r)) fail(func, "値が型に収まらない（答えを mod で求めるなら math_v2 の ModInt998 などを使う）");
        return r;
    } else {
        return a + b;
    }
}
// 数え上げの途中の足し算：T が整数型なら溢れたことを印 oa に残すだけ（途中で溢れても答えに使わなければ止めない。最後に印を見て止める）
template <class T>
void add_count(T& a, bool& oa, const T& b, bool ob) {
    if constexpr (std::is_integral<T>::value) {
        T r;
        oa = __builtin_add_overflow(a, b, &r) || oa || ob;
        a = r;
    } else {
        a = a + b;
    }
}
[[noreturn]] inline void count_overflow(const char* func) { fail(func, "値が型に収まらない（答えを mod で求めるなら math_v2 の ModInt998 などを使う）"); }
// 列 → 値の列（char は 0〜255）
template <class Seq>
std::vector<long long> codes_of(const Seq& s0) {
    const auto& s = seq(s0);
    std::vector<long long> c(s.size());
    for (std::size_t i = 0; i < s.size(); ++i) c[i] = code_of(s[i]);
    return c;
}
// 辞書順の比較（-1 / 0 / 1）。char は unsigned char として比べる（std::string と同じ順）
template <class A, class B>
int lex_compare(const A& a, const B& b) {
    const std::size_t n = a.size(), m = b.size(), k = n < m ? n : m;
    for (std::size_t i = 0; i < k; ++i) {
        const long long x = code_of(a[i]), y = code_of(b[i]);
        if (x != y) return x < y ? -1 : 1;
    }
    return n == m ? 0 : (n < m ? -1 : 1);
}
}  // namespace str_internal

// ---- 部分文字列の位置 [l, r) ----
struct StrSpan {
    int l = 0, r = 0;
    int len() const { return r - l; }
    bool empty() const { return r <= l; }
    friend bool operator==(const StrSpan& a, const StrSpan& b) { return a.l == b.l && a.r == b.r; }
    friend bool operator!=(const StrSpan& a, const StrSpan& b) { return !(a == b); }
    friend std::ostream& operator<<(std::ostream& os, const StrSpan& sp) { return os << "[" << sp.l << ", " << sp.r << ")"; }  // デバッグ用
};
template <class Seq>
str_internal::seq_t<Seq> substr_of(const Seq& s0, StrSpan sp) {  // s[sp.l, sp.r) のコピー（string なら string、vector なら vector）
    const auto& s = str_internal::seq(s0);
    str_internal::check_range("substr_of", sp.l, sp.r, (long long)s.size());
    return str_internal::seq_t<Seq>(s.begin() + sp.l, s.begin() + sp.r);
}

// ---- 文字の種類（Alphabet）：文字・値 → 番号 0〜size()-1。番号の順が辞書順として扱われる ----
//   Alphabet::lower()（a〜z。既定）/ upper()（A〜Z）/ digits()（0〜9）/ letters()（A〜Z, a〜z の 52 種。ASCII 順）
//   alnum()（0〜9, A〜Z, a〜z の 62 種。ASCII 順）/ ascii()（0〜127）/ bytes()（0〜255。日本語など 128 以上の char も）
//   of("()[]")（並べた文字だけ。並べた順が番号の順 = 辞書順）/ ints(K)（整数 0〜K-1）/ ints(lo, hi)（整数 [lo, hi)）
//   範囲外の文字・値を id に渡すとメッセージを出して止まる（contains で先に確かめられる）
struct Alphabet {
    Alphabet() : Alphabet(of("abcdefghijklmnopqrstuvwxyz")) {}
    static Alphabet lower() { return of("abcdefghijklmnopqrstuvwxyz"); }
    static Alphabet upper() { return of("ABCDEFGHIJKLMNOPQRSTUVWXYZ"); }
    static Alphabet digits() { return of("0123456789"); }
    static Alphabet letters() { return of("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"); }
    static Alphabet alnum() { return of("0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"); }
    static Alphabet ascii() { return code_range(0, 128); }
    static Alphabet bytes() { return code_range(0, 256); }
    static Alphabet of(const std::string& cs) {  // 並べた文字（同じ文字を 2 回並べると止まる）
        if (cs.empty()) str_internal::fail("Alphabet::of", "文字が 1 つもない");
        Alphabet A(0);
        A.kind_ = CHARS;
        for (std::size_t i = 0; i < cs.size(); ++i) {
            const unsigned char c = (unsigned char)cs[i];
            if (A.table_[c] >= 0) str_internal::fail("Alphabet::of", std::string("同じ文字が 2 回ある: '") + cs[i] + "'");
            A.table_[c] = (short)i;
            A.chars_[i] = c;
        }
        A.size_ = (int)cs.size();
        return A;
    }
    static Alphabet ints(long long K) { return ints(0, K); }  // 整数 0〜K-1
    static Alphabet ints(long long lo, long long hi) {         // 整数 [lo, hi)
        if (hi <= lo || (unsigned long long)hi - (unsigned long long)lo > 1000000000ULL)
            str_internal::fail("Alphabet::ints", "lo < hi かつ hi - lo ≤ 10^9 で指定する（大きな値は compress_values で 0〜K-1 に）");
        Alphabet A(0);
        A.kind_ = INTS;
        A.lo_ = lo;
        A.size_ = (int)(hi - lo);
        return A;
    }

    int size() const { return size_; }
    bool is_chars() const { return kind_ == CHARS; }  // 文字の表で決めたものか（of・lower など）。false なら ints
    template <class T>
    bool contains(T x) const { return index_of(str_internal::code_of(x)) >= 0; }
    template <class T>
    int id(T x) const {  // 文字・値の番号（範囲外は止める）
        const long long v = str_internal::code_of(x);
        const int r = index_of(v);
        if (r < 0) out_of_range(v);
        return r;
    }
    long long value(int i) const {  // 番号 i の値（文字なら文字コード 0〜255）
        check_index(i);
        return kind_ == CHARS ? (long long)chars_[i] : lo_ + i;
    }
    char ch(int i) const {  // 番号 i の文字（ints で 0〜255 の外の値なら止める）
        const long long v = value(i);
        if (v < 0 || v > 255) str_internal::fail("Alphabet::ch", "番号 " + std::to_string(i) + " の値 " + std::to_string(v) + " は文字にできない");
        return (char)(unsigned char)v;
    }
    std::string chars() const {  // 文字を番号の順に並べたもの（ints なら 0〜255 の範囲のときだけ）
        std::string s;
        for (int i = 0; i < size_; ++i) s += ch(i);
        return s;
    }

   private:
    enum Kind : unsigned char { CHARS, INTS };
    explicit Alphabet(int) { table_.fill(-1); }
    static Alphabet code_range(int lo, int hi) {
        Alphabet A(0);
        A.kind_ = CHARS;
        for (int c = lo; c < hi; ++c) A.table_[c] = (short)(c - lo), A.chars_[c - lo] = (unsigned char)c;
        A.size_ = hi - lo;
        return A;
    }
    int index_of(long long v) const {
        if (kind_ == CHARS) return (unsigned long long)v < 256ULL ? table_[(std::size_t)v] : -1;
        const unsigned long long d = (unsigned long long)v - (unsigned long long)lo_;
        return d < (unsigned long long)size_ ? (int)d : -1;
    }
    void check_index(int i) const {
        if (i < 0 || i >= size_) str_internal::fail("Alphabet", "番号 " + std::to_string(i) + " は 0〜" + std::to_string(size_ - 1) + " の外");
    }
    [[noreturn]] void out_of_range(long long v) const {
        std::string what = "値 " + std::to_string(v) + " は";
        if (kind_ == CHARS && v >= 32 && v < 127) what = std::string("文字 '") + (char)v + "'（コード " + std::to_string(v) + "）は";
        str_internal::fail("Alphabet::id", what + "文字の種類に入っていない（種類の数 " + std::to_string(size_) + "）");
    }
    Kind kind_ = CHARS;
    int size_ = 0;
    long long lo_ = 0;
    std::array<short, 256> table_{};
    std::array<unsigned char, 256> chars_{};
};

template <class Seq>
std::vector<int> to_ids(const Seq& s0, const Alphabet& A = Alphabet::lower()) {  // 各文字の番号の列
    const auto& s = str_internal::seq(s0);
    std::vector<int> r(s.size());
    for (std::size_t i = 0; i < s.size(); ++i) r[i] = A.id(s[i]);
    return r;
}
inline std::string from_ids(const std::vector<int>& ids, const Alphabet& A = Alphabet::lower()) {  // 番号の列 → 文字列
    std::string s(ids.size(), '\0');
    for (std::size_t i = 0; i < ids.size(); ++i) s[i] = A.ch(ids[i]);
    return s;
}

// ---- 数の列を 0〜K-1 に圧縮（大小の順を保つ）：Alphabet::ints(K) と組み合わせる ----
struct CompressedValues {
    std::vector<int> ids;            // 各要素の番号（0〜K-1。値の小さい順）
    std::vector<long long> values;   // 番号 → 元の値（昇順。K = values.size()）
    int size() const { return (int)values.size(); }
    Alphabet alphabet() const { return Alphabet::ints(values.empty() ? 1 : (long long)values.size()); }
};
template <class Seq>
CompressedValues compress_values(const Seq& s0) {
    const auto& s = str_internal::seq(s0);
    const std::size_t n = s.size();
    std::vector<std::pair<long long, int>> vi(n);
    for (std::size_t i = 0; i < n; ++i) vi[i] = {str_internal::code_of(s[i]), (int)i};
    std::sort(vi.begin(), vi.end());
    CompressedValues c;
    c.ids.resize(n);
    for (std::size_t i = 0; i < n; ++i) {
        if (i == 0 || vi[i].first != vi[i - 1].first) c.values.push_back(vi[i].first);
        c.ids[vi[i].second] = (int)c.values.size() - 1;
    }
    return c;
}

// ---- 文字の数 ----
template <class Seq>
std::vector<int> count_chars(const Seq& s0, const Alphabet& A = Alphabet::lower()) {  // 番号ごとの個数（長さ A.size()）
    const auto& s = str_internal::seq(s0);
    std::vector<int> c(A.size(), 0);
    for (std::size_t i = 0; i < s.size(); ++i) ++c[A.id(s[i])];
    return c;
}

// ---- 区間の文字の数・次 / 前の出現を O(1) で（前処理 O(Nσ)・メモリ (N+1)σ 個の int。N = 2×10^5, σ = 26 で約 21MB） ----
//   メモリを抑えたいとき・σ が大きいときは CharPositions（O(log N)・メモリ O(N + σ)）
struct CharPrefix {
    CharPrefix() = default;
    template <class Seq>
    explicit CharPrefix(const Seq& s0, const Alphabet& A = Alphabet::lower()) : A_(A) {
        const auto& s = str_internal::seq(s0);
        n_ = (int)s.size(), sg_ = A_.size();
        if ((long long)(n_ + 1) * sg_ > 600000000LL) str_internal::fail("CharPrefix", "(N + 1)σ が 6×10^8 を超える（メモリが足りない。CharPositions を使う）");
        cnt_.assign((std::size_t)(n_ + 1) * sg_, 0);
        std::vector<int> id(n_);
        start_.assign(sg_ + 1, 0);
        for (int i = 0; i < n_; ++i) {
            const int* src = cnt_.data() + (std::size_t)i * sg_;
            int* dst = cnt_.data() + (std::size_t)(i + 1) * sg_;
            std::copy(src, src + sg_, dst);
            id[i] = A_.id(s[i]);
            ++dst[id[i]], ++start_[id[i] + 1];
        }
        for (int k = 0; k < sg_; ++k) start_[k + 1] += start_[k];
        pos_.resize(n_);
        std::vector<int> cur(start_.begin(), start_.end() - 1);
        for (int i = 0; i < n_; ++i) pos_[cur[id[i]]++] = i;
    }
    int size() const { return n_; }
    const Alphabet& alphabet() const { return A_; }
    template <class C>
    int count(int l, int r, C c) const {  // [l, r) にある c の数  O(1)
        str_internal::check_range("CharPrefix::count", l, r, n_);
        const int k = A_.id(c);
        return cnt_[(std::size_t)r * sg_ + k] - cnt_[(std::size_t)l * sg_ + k];
    }
    template <class C>
    bool has(int l, int r, C c) const { return count(l, r, c) > 0; }  // [l, r) に c があるか  O(1)
    std::vector<int> counts(int l, int r) const {  // [l, r) の番号ごとの数  O(σ)
        str_internal::check_range("CharPrefix::counts", l, r, n_);
        std::vector<int> c(sg_);
        for (int k = 0; k < sg_; ++k) c[k] = cnt_[(std::size_t)r * sg_ + k] - cnt_[(std::size_t)l * sg_ + k];
        return c;
    }
    int distinct(int l, int r) const {  // [l, r) にある文字の種類の数  O(σ)
        str_internal::check_range("CharPrefix::distinct", l, r, n_);
        int d = 0;
        for (int k = 0; k < sg_; ++k) d += cnt_[(std::size_t)r * sg_ + k] != cnt_[(std::size_t)l * sg_ + k];
        return d;
    }
    template <class C>
    int next(int i, C c) const {  // i 以上で最初の c の位置（なければ -1。0 ≤ i ≤ N）  O(1)
        check_pos("CharPrefix::next", i);
        const int k = A_.id(c), m = cnt_[(std::size_t)i * sg_ + k];
        return start_[k] + m < start_[k + 1] ? pos_[start_[k] + m] : -1;
    }
    template <class C>
    int prev(int i, C c) const {  // i より前（[0, i)）で最後の c の位置（なければ -1。0 ≤ i ≤ N）  O(1)
        check_pos("CharPrefix::prev", i);
        const int k = A_.id(c), m = cnt_[(std::size_t)i * sg_ + k];
        return m > 0 ? pos_[start_[k] + m - 1] : -1;
    }
    template <class C>
    int kth(C c, long long k) const {  // c の k 番目（0 から）の位置（なければ -1）  O(1)
        const int id = A_.id(c);
        return (k >= 0 && k < start_[id + 1] - start_[id]) ? pos_[start_[id] + k] : -1;
    }
    template <class C>
    int total(C c) const {  // 全体での c の数  O(1)
        const int k = A_.id(c);
        return start_[k + 1] - start_[k];
    }

   private:
    void check_pos(const char* func, int i) const {
        if (i < 0 || i > n_) str_internal::fail(func, "位置 " + std::to_string(i) + " は 0〜N（= " + std::to_string(n_) + "）の外");
    }
    Alphabet A_;
    int n_ = 0, sg_ = 0;
    std::vector<int> cnt_;           // cnt_[i * σ + k] = [0, i) にある番号 k の数
    std::vector<int> start_, pos_;   // 番号ごとの出現位置（pos_[start_[k] .. start_[k+1])）
};

// ---- 文字の出現位置（前処理 O(N + σ)・メモリ O(N + σ)。σ が大きくても使える） ----
struct CharPositions {
    CharPositions() = default;
    template <class Seq>
    explicit CharPositions(const Seq& s0, const Alphabet& A = Alphabet::lower()) : A_(A) {
        const auto& s = str_internal::seq(s0);
        n_ = (int)s.size();
        const int sg = A_.size();
        std::vector<int> id(n_);
        start_.assign(sg + 1, 0);
        for (int i = 0; i < n_; ++i) id[i] = A_.id(s[i]), ++start_[id[i] + 1];
        for (int k = 0; k < sg; ++k) start_[k + 1] += start_[k];
        pos_.resize(n_), occ_.resize(n_);
        std::vector<int> cur(start_.begin(), start_.end() - 1);
        for (int i = 0; i < n_; ++i) occ_[i] = cur[id[i]] - start_[id[i]], pos_[cur[id[i]]++] = i;
    }
    int size() const { return n_; }
    const Alphabet& alphabet() const { return A_; }
    template <class C>
    int next(int i, C c) const {  // i 以上で最初の c の位置（なければ -1。0 ≤ i ≤ N）  O(log N)
        check_pos("CharPositions::next", i);
        const int k = A_.id(c), j = str_internal::lower_index(pos_.data(), start_[k], start_[k + 1], i);
        return j < start_[k + 1] ? pos_[j] : -1;
    }
    template <class C>
    int prev(int i, C c) const {  // i より前（[0, i)）で最後の c の位置（なければ -1。0 ≤ i ≤ N）  O(log N)
        check_pos("CharPositions::prev", i);
        const int k = A_.id(c), j = str_internal::lower_index(pos_.data(), start_[k], start_[k + 1], i);
        return j > start_[k] ? pos_[j - 1] : -1;
    }
    template <class C>
    int kth(C c, long long k) const {  // c の k 番目（0 から）の位置（なければ -1）  O(1)
        const int id = A_.id(c);
        return (k >= 0 && k < start_[id + 1] - start_[id]) ? pos_[start_[id] + k] : -1;
    }
    template <class C>
    int count(int l, int r, C c) const {  // [l, r) にある c の数  O(log N)
        str_internal::check_range("CharPositions::count", l, r, n_);
        const int k = A_.id(c);
        return str_internal::lower_index(pos_.data(), start_[k], start_[k + 1], r) - str_internal::lower_index(pos_.data(), start_[k], start_[k + 1], l);
    }
    template <class C>
    bool has(int l, int r, C c) const { return count(l, r, c) > 0; }
    template <class C>
    int total(C c) const {  // 全体での c の数  O(1)
        const int k = A_.id(c);
        return start_[k + 1] - start_[k];
    }
    template <class C>
    std::vector<int> positions(C c) const {  // c の位置の一覧（昇順）
        const int k = A_.id(c);
        return std::vector<int>(pos_.begin() + start_[k], pos_.begin() + start_[k + 1]);
    }
    int occurrence_index(int i) const {  // s[i] が「その文字の何番目（0 から）の出現」か  O(1)
        if (i < 0 || i >= n_) str_internal::fail("CharPositions::occurrence_index", "位置 " + std::to_string(i) + " は 0〜N-1 の外");
        return occ_[i];
    }

   private:
    void check_pos(const char* func, int i) const {
        if (i < 0 || i > n_) str_internal::fail(func, "位置 " + std::to_string(i) + " は 0〜N（= " + std::to_string(n_) + "）の外");
    }
    Alphabet A_;
    int n_ = 0;
    std::vector<int> start_, pos_, occ_;
};

// ---- ランレングス圧縮：連続する同じ値を {value, len, start} にまとめる ----
template <class T>
struct RunBlock {
    T value{};
    long long len = 0;    // 続く長さ
    long long start = 0;  // 元の列での開始位置
    friend std::ostream& operator<<(std::ostream& os, const RunBlock& b) { return os << "{" << b.value << ", " << b.len << ", " << b.start << "}"; }  // デバッグ用
};
template <class Seq>
std::vector<RunBlock<str_internal::elem_t<Seq>>> run_length(const Seq& s0) {
    const auto& s = str_internal::seq(s0);
    std::vector<RunBlock<str_internal::elem_t<Seq>>> r;
    std::size_t cnt = s.empty() ? 0 : 1;
    for (std::size_t i = 1; i < s.size(); ++i) cnt += !(s[i] == s[i - 1]);
    r.reserve(cnt);
    for (std::size_t i = 0, n = s.size(); i < n;) {
        std::size_t j = i + 1;
        while (j < n && s[j] == s[i]) ++j;
        r.push_back({s[i], (long long)(j - i), (long long)i});
        i = j;
    }
    return r;
}
template <class T>
auto run_length_decode(const std::vector<RunBlock<T>>& runs) {  // 元に戻す（T = char なら string、それ以外は vector<T>）
    using R = std::conditional_t<std::is_same<T, char>::value, std::string, std::vector<T>>;
    long long total = 0;
    for (const auto& b : runs) {
        if (b.len < 0) str_internal::fail("run_length_decode", "len が負");
        total += b.len;
        if (total > 2000000000LL) str_internal::fail("run_length_decode", "元の長さが 2×10^9 を超える");
    }
    R r;
    r.reserve((std::size_t)total);
    for (const auto& b : runs) r.insert(r.end(), (std::size_t)b.len, b.value);
    return r;
}

// ---- 変換 ----
inline std::string to_lower(std::string s) {  // 英大文字を小文字に（ほかはそのまま）
    for (char& c : s)
        if ('A' <= c && c <= 'Z') c = (char)(c - 'A' + 'a');
    return s;
}
inline std::string to_upper(std::string s) {
    for (char& c : s)
        if ('a' <= c && c <= 'z') c = (char)(c - 'a' + 'A');
    return s;
}
inline std::string swap_case(std::string s) {  // 大文字と小文字を入れ替える
    for (char& c : s) {
        if ('A' <= c && c <= 'Z') c = (char)(c - 'A' + 'a');
        else if ('a' <= c && c <= 'z') c = (char)(c - 'a' + 'A');
    }
    return s;
}
template <class Seq>
str_internal::seq_t<Seq> caesar_shift(const Seq& s0, long long k, const Alphabet& A = Alphabet::lower()) {  // A の文字を番号で k ずらす（負も。A にない文字はそのまま）
    auto s = str_internal::seq_t<Seq>(str_internal::seq(s0));
    const long long sg = A.size(), sh = ((k % sg) + sg) % sg;
    for (auto& c : s)
        if (A.contains(c)) c = (std::decay_t<decltype(c)>)A.value((int)((A.id(c) + sh) % sg));
    return s;
}
// 文字の置き換えを重ねる：replace(c, d) で「今 c になっている文字を全部 d に」。apply(s) で s に当てはめる（ABC342 C）
struct CharMap {
    explicit CharMap(const Alphabet& A = Alphabet::lower()) : A_(A), f_(A.size()) {
        for (int i = 0; i < A_.size(); ++i) f_[i] = i;
    }
    template <class C, class D>
    void replace(C c, D d) {  // 今 c になっているものを全部 d に  O(σ)
        const int x = A_.id(c), y = A_.id(d);
        if (x == y) return;
        for (int& v : f_)
            if (v == x) v = y;
    }
    template <class C, class D>
    void set(C c, D d) { f_[A_.id(c)] = A_.id(d); }  // 元の文字 c の行き先を d にする（ほかは変えない）
    template <class C>
    auto operator()(C c) const {  // 元の文字 c の行き先（A にない文字はそのまま）
        return A_.contains(c) ? (C)A_.value(f_[A_.id(c)]) : c;
    }
    template <class Seq>
    str_internal::seq_t<Seq> apply(const Seq& s0) const {  // s の各文字を行き先に（A にない文字はそのまま）
        auto s = str_internal::seq_t<Seq>(str_internal::seq(s0));
        for (auto& c : s) c = (*this)(c);
        return s;
    }
    CharMap compose(const CharMap& g) const {  // 「この置き換えのあとに g」を 1 つにしたもの
        CharMap h(A_);
        for (int i = 0; i < A_.size(); ++i) h.f_[i] = g.f_[f_[i]];
        return h;
    }
    const std::vector<int>& table() const { return f_; }  // 番号 → 行き先の番号

   private:
    Alphabet A_;
    std::vector<int> f_;
};
template <class Seq>
str_internal::seq_t<Seq> reversed(const Seq& s0) {  // 反転したコピー
    const auto& s = str_internal::seq(s0);
    return str_internal::seq_t<Seq>(s.rbegin(), s.rend());
}
template <class Seq>
str_internal::seq_t<Seq> min_with_reverse(const Seq& s0) {  // s と反転のうち辞書順で小さい方（向きを区別しないときの代表。ABC310 C）
    auto s = str_internal::seq_t<Seq>(str_internal::seq(s0));
    auto t = reversed(s);
    return str_internal::lex_compare(t, s) < 0 ? t : s;
}
template <class Seq>
str_internal::seq_t<Seq> sorted_chars(const Seq& s0) {  // 並べ替えたもの（並べ替えて同じかの比較用）
    auto s = str_internal::seq_t<Seq>(str_internal::seq(s0));
    std::sort(s.begin(), s.end(), [](const auto& a, const auto& b) { return str_internal::code_of(a) < str_internal::code_of(b); });
    return s;
}
template <class Seq>
str_internal::seq_t<Seq> repeat_string(const Seq& s0, long long k) {  // s を k 回つなげる
    const auto& s = str_internal::seq(s0);
    if (k < 0) str_internal::fail("repeat_string", "k ≥ 0 で使う");
    if (!s.empty() && k > 2000000000LL / (long long)s.size()) str_internal::fail("repeat_string", "長さが 2×10^9 を超える");
    str_internal::seq_t<Seq> r;
    r.reserve((std::size_t)(k * (long long)s.size()));
    for (long long i = 0; i < k; ++i) r.insert(r.end(), s.begin(), s.end());
    return r;
}

// ---- 比較 ----
template <class A, class B>
int common_prefix_length(const A& a0, const B& b0) {  // 共通の接頭辞の長さ
    const auto& a = str_internal::seq(a0);
    const auto& b = str_internal::seq(b0);
    const std::size_t k = a.size() < b.size() ? a.size() : b.size();
    std::size_t i = 0;
    while (i < k && str_internal::code_of(a[i]) == str_internal::code_of(b[i])) ++i;
    return (int)i;
}
template <class A, class B>
int common_suffix_length(const A& a0, const B& b0) {  // 共通の接尾辞の長さ
    const auto& a = str_internal::seq(a0);
    const auto& b = str_internal::seq(b0);
    const std::size_t k = a.size() < b.size() ? a.size() : b.size();
    std::size_t i = 0;
    while (i < k && str_internal::code_of(a[a.size() - 1 - i]) == str_internal::code_of(b[b.size() - 1 - i])) ++i;
    return (int)i;
}
template <class A, class B>
int hamming_distance(const A& a0, const B& b0) {  // 同じ位置で違う文字の数（長さが同じこと）
    const auto& a = str_internal::seq(a0);
    const auto& b = str_internal::seq(b0);
    if (a.size() != b.size()) str_internal::fail("hamming_distance", "長さが違う");
    int d = 0;
    for (std::size_t i = 0; i < a.size(); ++i) d += str_internal::code_of(a[i]) != str_internal::code_of(b[i]);
    return d;
}
template <class T, class S>
bool is_subsequence(const T& t0, const S& s0) {  // t が s の部分列（飛び飛びでよい）か  O(|s|)。何度も聞くなら 06_subsequence の NextTable
    const auto& t = str_internal::seq(t0);
    const auto& s = str_internal::seq(s0);
    std::size_t j = 0;
    for (std::size_t i = 0; i < s.size() && j < t.size(); ++i)
        if (str_internal::code_of(s[i]) == str_internal::code_of(t[j])) ++j;
    return j == t.size();
}
template <class Seq>
bool is_palindrome(const Seq& s0) {  // 全体が回文か  O(N)。区間を何度も聞くなら 04_palindrome の Palindromes
    const auto& s = str_internal::seq(s0);
    for (std::size_t i = 0, n = s.size(); i < n / 2; ++i)
        if (!(s[i] == s[n - 1 - i])) return false;
    return true;
}

// ---- 区切り・つなげる ----
inline std::vector<std::string> split(const std::string& s, char sep, bool skip_empty = false) {  // sep で区切る（空の部分も残す。skip_empty で除く）
    std::vector<std::string> r;
    std::size_t b = 0;
    for (std::size_t i = 0; i <= s.size(); ++i)
        if (i == s.size() || s[i] == sep) {
            if (!(skip_empty && i == b)) r.push_back(s.substr(b, i - b));
            b = i + 1;
        }
    return r;
}
template <class V>
std::string join(const std::vector<V>& vs, const std::string& sep = " ") {  // 区切りを挟んでつなげる（数も文字列にする。答えの出力に）
    std::string r;
    for (std::size_t i = 0; i < vs.size(); ++i) {
        if (i) r += sep;
        if constexpr (std::is_same<V, std::string>::value) r += vs[i];
        else if constexpr (std::is_same<V, char>::value) r += vs[i];
        else if constexpr (std::is_integral<V>::value) r += std::to_string(vs[i]);
        else {
            std::ostringstream os;
            os << vs[i];
            r += os.str();
        }
    }
    return r;
}
#endif  // GG_STR_BASE
#ifndef GG_LIS
#define GG_LIS
// ---- 最長増加部分列（LIS。文字列第二世代） ----
// 列は vector<整数>（int・long long。負・10^18 級もよい）でも string でもよい。strict = true で狭義（<）、false で広義（≤）
// 位置・添字は 0 始まり。「部分列」は飛び飛びに選んだもの（順番はそのまま）
namespace str_internal {
// a[0, len) は単調増加。strict なら a[k] ≥ x の最初の k、そうでなければ a[k] > x の最初の k
// （末尾より大きければすぐ返す。ばらばらな列では分岐ありの二分探索の方が速かった）
inline int lis_slot(const long long* a, int len, long long x, bool strict) {
    if (len == 0 || (strict ? a[len - 1] < x : a[len - 1] <= x)) return len;
    int lo = -1, hi = len - 1;
    if (strict) {
        while (hi - lo > 1) {
            const int mid = (lo + hi) >> 1;
            if (a[mid] < x) lo = mid;
            else hi = mid;
        }
    } else {
        while (hi - lo > 1) {
            const int mid = (lo + hi) >> 1;
            if (a[mid] <= x) lo = mid;
            else hi = mid;
        }
    }
    return hi;
}
// end[i] = i で終わる最長の長さ。prev を渡すと復元用の 1 つ前の添字も、last には最長の列の最後の添字
inline int lis_core(const std::vector<long long>& v, bool strict, std::vector<int>* end, std::vector<int>* prev, int* last = nullptr) {
    const int n = (int)v.size();
    std::vector<long long> tail(n);
    std::vector<int> tail_id(prev || last ? n : 0);
    if (end) end->resize(n);
    if (prev) prev->resize(n);
    int L = 0;
    for (int i = 0; i < n; ++i) {
        const int p = lis_slot(tail.data(), L, v[i], strict);
        tail[p] = v[i];
        if (end) (*end)[i] = p + 1;
        if (prev) (*prev)[i] = p > 0 ? tail_id[p - 1] : -1;
        if (prev || last) tail_id[p] = i;
        if (p == L) ++L;
    }
    if (last) *last = L > 0 ? tail_id[L - 1] : -1;
    return L;
}
// 大小を逆にする（~v = -v - 1 は溢れない）
inline std::vector<long long> flipped(std::vector<long long> v) {
    for (auto& x : v) x = ~x;
    return v;
}
}  // namespace str_internal

template <class Seq>
int lis_length(const Seq& a, bool strict = true) {  // 最長増加部分列の長さ  O(N log N)
    return str_internal::lis_core(str_internal::codes_of(a), strict, nullptr, nullptr);
}
template <class Seq>
std::vector<int> lis_restore(const Seq& a, bool strict = true) {  // 最長増加部分列を 1 つ（添字の列。昇順）  O(N log N)
    std::vector<int> prev;
    int cur = -1;
    const int L = str_internal::lis_core(str_internal::codes_of(a), strict, nullptr, &prev, &cur);
    std::vector<int> r(L);
    for (int k = L - 1; k >= 0; --k) r[k] = cur, cur = prev[cur];
    return r;
}
template <class Seq>
std::vector<int> lis_ending_at(const Seq& a, bool strict = true) {  // r[i] = a[i] で終わる最長の増加部分列の長さ  O(N log N)
    std::vector<int> end;
    str_internal::lis_core(str_internal::codes_of(a), strict, &end, nullptr);
    return end;
}
template <class Seq>
std::vector<int> lis_starting_at(const Seq& a, bool strict = true) {  // r[i] = a[i] から始まる最長の増加部分列の長さ  O(N log N)
    auto v = str_internal::flipped(str_internal::codes_of(a));
    std::reverse(v.begin(), v.end());  // 後ろから読むと「減少」なので大小を逆にして増加に
    std::vector<int> end;
    str_internal::lis_core(v, strict, &end, nullptr);
    std::reverse(end.begin(), end.end());
    return end;
}

struct LisMembers {
    int length = 0;                // LIS の長さ
    std::vector<bool> in_some;     // in_some[i] = a[i] を含む LIS がある（ABC354 F）
    std::vector<bool> in_all;      // in_all[i]  = すべての LIS が a[i] を含む
};
template <class Seq>
LisMembers lis_members(const Seq& a, bool strict = true) {  // 各位置が LIS に入りうるか・必ず入るか  O(N log N)
    const auto e = lis_ending_at(a, strict), s = lis_starting_at(a, strict);
    const int n = (int)e.size();
    LisMembers r;
    for (int i = 0; i < n; ++i) r.length = std::max(r.length, e[i]);
    r.in_some.assign(n, false), r.in_all.assign(n, false);
    std::vector<int> cnt(r.length + 1, 0);  // 長さの段 k（LIS の k 番目）に入りうる位置の数
    for (int i = 0; i < n; ++i)
        if (e[i] + s[i] - 1 == r.length) r.in_some[i] = true, ++cnt[e[i]];
    for (int i = 0; i < n; ++i) r.in_all[i] = r.in_some[i] && cnt[e[i]] == 1;
    return r;
}

template <class T, class Seq>
T count_lis(const Seq& a, bool strict = true) {  // LIS の数（添字の列として数える。T は ModInt998・long long（答えが溢れたら止める）・double など）  O(N log N)
    const CompressedValues cv = compress_values(a);
    const int n = (int)cv.ids.size(), K = cv.size();
    std::vector<int> len(K + 1, 0);  // フェニック木：値の番号の接頭辞での（最長の長さ, その数, 溢れたか）
    std::vector<T> cnt(K + 1, T(0));
    std::vector<char> over(K + 1, 0);
    auto merge = [&](int& l0, T& c0, bool& o0, int l1, const T& c1, bool o1) {
        if (l1 > l0) l0 = l1, c0 = c1, o0 = o1;
        else if (l1 == l0 && l1 > 0) str_internal::add_count(c0, o0, c1, o1);
    };
    for (int i = 0; i < n; ++i) {
        const int id = cv.ids[i];
        int bl = 0;
        T bc = T(0);
        bool bo = false;
        for (int k = strict ? id : id + 1; k > 0; k -= k & -k) merge(bl, bc, bo, len[k], cnt[k], over[k]);
        if (bl == 0) bc = T(1);
        ++bl;
        for (int k = id + 1; k <= K; k += k & -k) {
            bool o = over[k];
            merge(len[k], cnt[k], o, bl, bc, bo);
            over[k] = o;
        }
    }
    int bl = 0;
    T bc = T(0);
    bool bo = false;
    for (int k = K; k > 0; k -= k & -k) merge(bl, bc, bo, len[k], cnt[k], over[k]);
    if (bo) str_internal::count_overflow("count_lis");
    return n == 0 ? T(1) : bc;  // 空の列の LIS は空の 1 通り
}

template <class Seq>
int lds_length(const Seq& a, bool strict = true) {  // 最長減少部分列の長さ  O(N log N)
    return str_internal::lis_core(str_internal::flipped(str_internal::codes_of(a)), strict, nullptr, nullptr);
}
template <class Seq>
int min_increasing_cover(const Seq& a, bool strict = true) {  // 増加部分列（strict なら狭義）に分ける最小の数 = 最長の「広義 / 狭義」減少部分列（Dilworth。ABC134 E）
    return lds_length(a, !strict);
}
template <class Seq>
int min_decreasing_cover(const Seq& a, bool strict = true) {  // 減少部分列に分ける最小の数
    return lis_length(a, !strict);
}

template <class X, class Y>
std::vector<int> lis_2d(const std::vector<std::pair<X, Y>>& ps) {  // x も y も狭義に増える最長の列（添字。x の小さい順。ABC038 D）  O(N log N)
    const int n = (int)ps.size();
    std::vector<int> ord(n);
    for (int i = 0; i < n; ++i) ord[i] = i;
    std::sort(ord.begin(), ord.end(), [&](int i, int j) {  // x の昇順、同じ x は y の降順（同じ x を 2 つ選ばないように）
        if (ps[i].first < ps[j].first) return true;
        if (ps[j].first < ps[i].first) return false;
        if (ps[j].second < ps[i].second) return true;
        if (ps[i].second < ps[j].second) return false;
        return i < j;
    });
    // y を < だけで順位に（double・unsigned long long などでも大小を保つ）
    std::vector<int> by(n), ry(n);
    for (int i = 0; i < n; ++i) by[i] = i;
    std::sort(by.begin(), by.end(), [&](int i, int j) { return ps[i].second < ps[j].second; });
    for (int k = 0, r = -1; k < n; ++k) {
        if (k == 0 || ps[by[k - 1]].second < ps[by[k]].second) ++r;
        ry[by[k]] = r;
    }
    std::vector<long long> ys(n);
    for (int k = 0; k < n; ++k) ys[k] = ry[ord[k]];
    std::vector<int> r = lis_restore(ys, true);
    for (int& k : r) k = ord[k];
    return r;
}
#endif  // GG_LIS

#ifndef GG_STR_BASE
#define GG_STR_BASE
// ---- 文字列基盤（文字列第二世代）  文字の種類 Alphabet・文字の数と位置・ランレングス・変換・比較・区切り ----
// 列（string / vector<整数>）を受け取る関数はどちらでもよい。文字列リテラル（"abc"）も渡せる
// 範囲は半開区間 [l, r)（0 始まり）。位置が「ない」ときは -1。使える範囲の外の引数はメッセージを出して止める
namespace str_internal {
[[noreturn]] inline void fail(const char* func, const std::string& msg) {
    std::cerr << func << ": " << msg << "\n";
    std::abort();
}
// 文字列リテラルを std::string に（それ以外はそのまま参照を返す）
template <class T>
const T& seq(const T& s) { return s; }
template <std::size_t N>
std::string seq(const char (&s)[N]) { return std::string(s); }
inline std::string seq(const char* s) { return std::string(s); }
template <class T>
using seq_t = std::decay_t<decltype(seq(std::declval<const T&>()))>;
template <class T>
using elem_t = std::decay_t<decltype(std::declval<const seq_t<T>&>()[0])>;
// 要素を「値」に（char は unsigned char として 0〜255。負の添字にしない）
template <class T>
long long code_of(T x) {
    if constexpr (std::is_same<T, char>::value || std::is_same<T, signed char>::value || std::is_same<T, unsigned char>::value)
        return (long long)(unsigned char)x;
    else
        return (long long)x;
}
// 範囲の確認（止めるときの文字列作りは別の関数に分けて、確認そのものを軽くする）
[[noreturn]] __attribute__((noinline, cold)) inline void range_fail(const char* func, long long l, long long r, long long n) {
    fail(func, "範囲 [l, r) は 0 ≤ l ≤ r ≤ 長さ で指定する（l = " + std::to_string(l) + ", r = " + std::to_string(r) + ", 長さ = " + std::to_string(n) + "）");
}
inline void check_range(const char* func, long long l, long long r, long long n) {
    if (__builtin_expect(l < 0 || r > n || l > r, 0)) range_fail(func, l, r, n);
}
// a[lo, hi) の中で x 以上の最初の位置（a は昇順）。分岐の少ない二分探索
inline int lower_index(const int* a, int lo, int hi, int x) {
    int len = hi - lo;
    if (len <= 0) return lo;
    const int* b = a + lo;
    while (len > 1) {
        const int half = len >> 1;
        b += (b[half - 1] < x) ? half : 0;
        len -= half;
    }
    return (int)(b - a) + (*b < x);
}
// 昇順の a の中で x より大きい最初の位置（std::upper_bound と同じ。デバッグビルドでも O(log N)）
template <class T>
int upper_index(const std::vector<T>& a, const T& x) {
    int lo = 0, hi = (int)a.size();
    while (lo < hi) {
        const int mid = (int)((unsigned)(lo + hi) >> 1);
        if (!(x < a[mid])) lo = mid + 1;
        else hi = mid;
    }
    return lo;
}
// 数え上げの足し算：T が整数型なら溢れたら止める（mint・double などはそのまま）
template <class T>
T checked_add(const T& a, const T& b, const char* func) {
    if constexpr (std::is_integral<T>::value) {
        T r;
        if (__builtin_add_overflow(a, b, &r)) fail(func, "値が型に収まらない（答えを mod で求めるなら math_v2 の ModInt998 などを使う）");
        return r;
    } else {
        return a + b;
    }
}
// 数え上げの途中の足し算：T が整数型なら溢れたことを印 oa に残すだけ（途中で溢れても答えに使わなければ止めない。最後に印を見て止める）
template <class T>
void add_count(T& a, bool& oa, const T& b, bool ob) {
    if constexpr (std::is_integral<T>::value) {
        T r;
        oa = __builtin_add_overflow(a, b, &r) || oa || ob;
        a = r;
    } else {
        a = a + b;
    }
}
[[noreturn]] inline void count_overflow(const char* func) { fail(func, "値が型に収まらない（答えを mod で求めるなら math_v2 の ModInt998 などを使う）"); }
// 列 → 値の列（char は 0〜255）
template <class Seq>
std::vector<long long> codes_of(const Seq& s0) {
    const auto& s = seq(s0);
    std::vector<long long> c(s.size());
    for (std::size_t i = 0; i < s.size(); ++i) c[i] = code_of(s[i]);
    return c;
}
// 辞書順の比較（-1 / 0 / 1）。char は unsigned char として比べる（std::string と同じ順）
template <class A, class B>
int lex_compare(const A& a, const B& b) {
    const std::size_t n = a.size(), m = b.size(), k = n < m ? n : m;
    for (std::size_t i = 0; i < k; ++i) {
        const long long x = code_of(a[i]), y = code_of(b[i]);
        if (x != y) return x < y ? -1 : 1;
    }
    return n == m ? 0 : (n < m ? -1 : 1);
}
}  // namespace str_internal

// ---- 部分文字列の位置 [l, r) ----
struct StrSpan {
    int l = 0, r = 0;
    int len() const { return r - l; }
    bool empty() const { return r <= l; }
    friend bool operator==(const StrSpan& a, const StrSpan& b) { return a.l == b.l && a.r == b.r; }
    friend bool operator!=(const StrSpan& a, const StrSpan& b) { return !(a == b); }
    friend std::ostream& operator<<(std::ostream& os, const StrSpan& sp) { return os << "[" << sp.l << ", " << sp.r << ")"; }  // デバッグ用
};
template <class Seq>
str_internal::seq_t<Seq> substr_of(const Seq& s0, StrSpan sp) {  // s[sp.l, sp.r) のコピー（string なら string、vector なら vector）
    const auto& s = str_internal::seq(s0);
    str_internal::check_range("substr_of", sp.l, sp.r, (long long)s.size());
    return str_internal::seq_t<Seq>(s.begin() + sp.l, s.begin() + sp.r);
}

// ---- 文字の種類（Alphabet）：文字・値 → 番号 0〜size()-1。番号の順が辞書順として扱われる ----
//   Alphabet::lower()（a〜z。既定）/ upper()（A〜Z）/ digits()（0〜9）/ letters()（A〜Z, a〜z の 52 種。ASCII 順）
//   alnum()（0〜9, A〜Z, a〜z の 62 種。ASCII 順）/ ascii()（0〜127）/ bytes()（0〜255。日本語など 128 以上の char も）
//   of("()[]")（並べた文字だけ。並べた順が番号の順 = 辞書順）/ ints(K)（整数 0〜K-1）/ ints(lo, hi)（整数 [lo, hi)）
//   範囲外の文字・値を id に渡すとメッセージを出して止まる（contains で先に確かめられる）
struct Alphabet {
    Alphabet() : Alphabet(of("abcdefghijklmnopqrstuvwxyz")) {}
    static Alphabet lower() { return of("abcdefghijklmnopqrstuvwxyz"); }
    static Alphabet upper() { return of("ABCDEFGHIJKLMNOPQRSTUVWXYZ"); }
    static Alphabet digits() { return of("0123456789"); }
    static Alphabet letters() { return of("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"); }
    static Alphabet alnum() { return of("0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"); }
    static Alphabet ascii() { return code_range(0, 128); }
    static Alphabet bytes() { return code_range(0, 256); }
    static Alphabet of(const std::string& cs) {  // 並べた文字（同じ文字を 2 回並べると止まる）
        if (cs.empty()) str_internal::fail("Alphabet::of", "文字が 1 つもない");
        Alphabet A(0);
        A.kind_ = CHARS;
        for (std::size_t i = 0; i < cs.size(); ++i) {
            const unsigned char c = (unsigned char)cs[i];
            if (A.table_[c] >= 0) str_internal::fail("Alphabet::of", std::string("同じ文字が 2 回ある: '") + cs[i] + "'");
            A.table_[c] = (short)i;
            A.chars_[i] = c;
        }
        A.size_ = (int)cs.size();
        return A;
    }
    static Alphabet ints(long long K) { return ints(0, K); }  // 整数 0〜K-1
    static Alphabet ints(long long lo, long long hi) {         // 整数 [lo, hi)
        if (hi <= lo || (unsigned long long)hi - (unsigned long long)lo > 1000000000ULL)
            str_internal::fail("Alphabet::ints", "lo < hi かつ hi - lo ≤ 10^9 で指定する（大きな値は compress_values で 0〜K-1 に）");
        Alphabet A(0);
        A.kind_ = INTS;
        A.lo_ = lo;
        A.size_ = (int)(hi - lo);
        return A;
    }

    int size() const { return size_; }
    bool is_chars() const { return kind_ == CHARS; }  // 文字の表で決めたものか（of・lower など）。false なら ints
    template <class T>
    bool contains(T x) const { return index_of(str_internal::code_of(x)) >= 0; }
    template <class T>
    int id(T x) const {  // 文字・値の番号（範囲外は止める）
        const long long v = str_internal::code_of(x);
        const int r = index_of(v);
        if (r < 0) out_of_range(v);
        return r;
    }
    long long value(int i) const {  // 番号 i の値（文字なら文字コード 0〜255）
        check_index(i);
        return kind_ == CHARS ? (long long)chars_[i] : lo_ + i;
    }
    char ch(int i) const {  // 番号 i の文字（ints で 0〜255 の外の値なら止める）
        const long long v = value(i);
        if (v < 0 || v > 255) str_internal::fail("Alphabet::ch", "番号 " + std::to_string(i) + " の値 " + std::to_string(v) + " は文字にできない");
        return (char)(unsigned char)v;
    }
    std::string chars() const {  // 文字を番号の順に並べたもの（ints なら 0〜255 の範囲のときだけ）
        std::string s;
        for (int i = 0; i < size_; ++i) s += ch(i);
        return s;
    }

   private:
    enum Kind : unsigned char { CHARS, INTS };
    explicit Alphabet(int) { table_.fill(-1); }
    static Alphabet code_range(int lo, int hi) {
        Alphabet A(0);
        A.kind_ = CHARS;
        for (int c = lo; c < hi; ++c) A.table_[c] = (short)(c - lo), A.chars_[c - lo] = (unsigned char)c;
        A.size_ = hi - lo;
        return A;
    }
    int index_of(long long v) const {
        if (kind_ == CHARS) return (unsigned long long)v < 256ULL ? table_[(std::size_t)v] : -1;
        const unsigned long long d = (unsigned long long)v - (unsigned long long)lo_;
        return d < (unsigned long long)size_ ? (int)d : -1;
    }
    void check_index(int i) const {
        if (i < 0 || i >= size_) str_internal::fail("Alphabet", "番号 " + std::to_string(i) + " は 0〜" + std::to_string(size_ - 1) + " の外");
    }
    [[noreturn]] void out_of_range(long long v) const {
        std::string what = "値 " + std::to_string(v) + " は";
        if (kind_ == CHARS && v >= 32 && v < 127) what = std::string("文字 '") + (char)v + "'（コード " + std::to_string(v) + "）は";
        str_internal::fail("Alphabet::id", what + "文字の種類に入っていない（種類の数 " + std::to_string(size_) + "）");
    }
    Kind kind_ = CHARS;
    int size_ = 0;
    long long lo_ = 0;
    std::array<short, 256> table_{};
    std::array<unsigned char, 256> chars_{};
};

template <class Seq>
std::vector<int> to_ids(const Seq& s0, const Alphabet& A = Alphabet::lower()) {  // 各文字の番号の列
    const auto& s = str_internal::seq(s0);
    std::vector<int> r(s.size());
    for (std::size_t i = 0; i < s.size(); ++i) r[i] = A.id(s[i]);
    return r;
}
inline std::string from_ids(const std::vector<int>& ids, const Alphabet& A = Alphabet::lower()) {  // 番号の列 → 文字列
    std::string s(ids.size(), '\0');
    for (std::size_t i = 0; i < ids.size(); ++i) s[i] = A.ch(ids[i]);
    return s;
}

// ---- 数の列を 0〜K-1 に圧縮（大小の順を保つ）：Alphabet::ints(K) と組み合わせる ----
struct CompressedValues {
    std::vector<int> ids;            // 各要素の番号（0〜K-1。値の小さい順）
    std::vector<long long> values;   // 番号 → 元の値（昇順。K = values.size()）
    int size() const { return (int)values.size(); }
    Alphabet alphabet() const { return Alphabet::ints(values.empty() ? 1 : (long long)values.size()); }
};
template <class Seq>
CompressedValues compress_values(const Seq& s0) {
    const auto& s = str_internal::seq(s0);
    const std::size_t n = s.size();
    std::vector<std::pair<long long, int>> vi(n);
    for (std::size_t i = 0; i < n; ++i) vi[i] = {str_internal::code_of(s[i]), (int)i};
    std::sort(vi.begin(), vi.end());
    CompressedValues c;
    c.ids.resize(n);
    for (std::size_t i = 0; i < n; ++i) {
        if (i == 0 || vi[i].first != vi[i - 1].first) c.values.push_back(vi[i].first);
        c.ids[vi[i].second] = (int)c.values.size() - 1;
    }
    return c;
}

// ---- 文字の数 ----
template <class Seq>
std::vector<int> count_chars(const Seq& s0, const Alphabet& A = Alphabet::lower()) {  // 番号ごとの個数（長さ A.size()）
    const auto& s = str_internal::seq(s0);
    std::vector<int> c(A.size(), 0);
    for (std::size_t i = 0; i < s.size(); ++i) ++c[A.id(s[i])];
    return c;
}

// ---- 区間の文字の数・次 / 前の出現を O(1) で（前処理 O(Nσ)・メモリ (N+1)σ 個の int。N = 2×10^5, σ = 26 で約 21MB） ----
//   メモリを抑えたいとき・σ が大きいときは CharPositions（O(log N)・メモリ O(N + σ)）
struct CharPrefix {
    CharPrefix() = default;
    template <class Seq>
    explicit CharPrefix(const Seq& s0, const Alphabet& A = Alphabet::lower()) : A_(A) {
        const auto& s = str_internal::seq(s0);
        n_ = (int)s.size(), sg_ = A_.size();
        if ((long long)(n_ + 1) * sg_ > 600000000LL) str_internal::fail("CharPrefix", "(N + 1)σ が 6×10^8 を超える（メモリが足りない。CharPositions を使う）");
        cnt_.assign((std::size_t)(n_ + 1) * sg_, 0);
        std::vector<int> id(n_);
        start_.assign(sg_ + 1, 0);
        for (int i = 0; i < n_; ++i) {
            const int* src = cnt_.data() + (std::size_t)i * sg_;
            int* dst = cnt_.data() + (std::size_t)(i + 1) * sg_;
            std::copy(src, src + sg_, dst);
            id[i] = A_.id(s[i]);
            ++dst[id[i]], ++start_[id[i] + 1];
        }
        for (int k = 0; k < sg_; ++k) start_[k + 1] += start_[k];
        pos_.resize(n_);
        std::vector<int> cur(start_.begin(), start_.end() - 1);
        for (int i = 0; i < n_; ++i) pos_[cur[id[i]]++] = i;
    }
    int size() const { return n_; }
    const Alphabet& alphabet() const { return A_; }
    template <class C>
    int count(int l, int r, C c) const {  // [l, r) にある c の数  O(1)
        str_internal::check_range("CharPrefix::count", l, r, n_);
        const int k = A_.id(c);
        return cnt_[(std::size_t)r * sg_ + k] - cnt_[(std::size_t)l * sg_ + k];
    }
    template <class C>
    bool has(int l, int r, C c) const { return count(l, r, c) > 0; }  // [l, r) に c があるか  O(1)
    std::vector<int> counts(int l, int r) const {  // [l, r) の番号ごとの数  O(σ)
        str_internal::check_range("CharPrefix::counts", l, r, n_);
        std::vector<int> c(sg_);
        for (int k = 0; k < sg_; ++k) c[k] = cnt_[(std::size_t)r * sg_ + k] - cnt_[(std::size_t)l * sg_ + k];
        return c;
    }
    int distinct(int l, int r) const {  // [l, r) にある文字の種類の数  O(σ)
        str_internal::check_range("CharPrefix::distinct", l, r, n_);
        int d = 0;
        for (int k = 0; k < sg_; ++k) d += cnt_[(std::size_t)r * sg_ + k] != cnt_[(std::size_t)l * sg_ + k];
        return d;
    }
    template <class C>
    int next(int i, C c) const {  // i 以上で最初の c の位置（なければ -1。0 ≤ i ≤ N）  O(1)
        check_pos("CharPrefix::next", i);
        const int k = A_.id(c), m = cnt_[(std::size_t)i * sg_ + k];
        return start_[k] + m < start_[k + 1] ? pos_[start_[k] + m] : -1;
    }
    template <class C>
    int prev(int i, C c) const {  // i より前（[0, i)）で最後の c の位置（なければ -1。0 ≤ i ≤ N）  O(1)
        check_pos("CharPrefix::prev", i);
        const int k = A_.id(c), m = cnt_[(std::size_t)i * sg_ + k];
        return m > 0 ? pos_[start_[k] + m - 1] : -1;
    }
    template <class C>
    int kth(C c, long long k) const {  // c の k 番目（0 から）の位置（なければ -1）  O(1)
        const int id = A_.id(c);
        return (k >= 0 && k < start_[id + 1] - start_[id]) ? pos_[start_[id] + k] : -1;
    }
    template <class C>
    int total(C c) const {  // 全体での c の数  O(1)
        const int k = A_.id(c);
        return start_[k + 1] - start_[k];
    }

   private:
    void check_pos(const char* func, int i) const {
        if (i < 0 || i > n_) str_internal::fail(func, "位置 " + std::to_string(i) + " は 0〜N（= " + std::to_string(n_) + "）の外");
    }
    Alphabet A_;
    int n_ = 0, sg_ = 0;
    std::vector<int> cnt_;           // cnt_[i * σ + k] = [0, i) にある番号 k の数
    std::vector<int> start_, pos_;   // 番号ごとの出現位置（pos_[start_[k] .. start_[k+1])）
};

// ---- 文字の出現位置（前処理 O(N + σ)・メモリ O(N + σ)。σ が大きくても使える） ----
struct CharPositions {
    CharPositions() = default;
    template <class Seq>
    explicit CharPositions(const Seq& s0, const Alphabet& A = Alphabet::lower()) : A_(A) {
        const auto& s = str_internal::seq(s0);
        n_ = (int)s.size();
        const int sg = A_.size();
        std::vector<int> id(n_);
        start_.assign(sg + 1, 0);
        for (int i = 0; i < n_; ++i) id[i] = A_.id(s[i]), ++start_[id[i] + 1];
        for (int k = 0; k < sg; ++k) start_[k + 1] += start_[k];
        pos_.resize(n_), occ_.resize(n_);
        std::vector<int> cur(start_.begin(), start_.end() - 1);
        for (int i = 0; i < n_; ++i) occ_[i] = cur[id[i]] - start_[id[i]], pos_[cur[id[i]]++] = i;
    }
    int size() const { return n_; }
    const Alphabet& alphabet() const { return A_; }
    template <class C>
    int next(int i, C c) const {  // i 以上で最初の c の位置（なければ -1。0 ≤ i ≤ N）  O(log N)
        check_pos("CharPositions::next", i);
        const int k = A_.id(c), j = str_internal::lower_index(pos_.data(), start_[k], start_[k + 1], i);
        return j < start_[k + 1] ? pos_[j] : -1;
    }
    template <class C>
    int prev(int i, C c) const {  // i より前（[0, i)）で最後の c の位置（なければ -1。0 ≤ i ≤ N）  O(log N)
        check_pos("CharPositions::prev", i);
        const int k = A_.id(c), j = str_internal::lower_index(pos_.data(), start_[k], start_[k + 1], i);
        return j > start_[k] ? pos_[j - 1] : -1;
    }
    template <class C>
    int kth(C c, long long k) const {  // c の k 番目（0 から）の位置（なければ -1）  O(1)
        const int id = A_.id(c);
        return (k >= 0 && k < start_[id + 1] - start_[id]) ? pos_[start_[id] + k] : -1;
    }
    template <class C>
    int count(int l, int r, C c) const {  // [l, r) にある c の数  O(log N)
        str_internal::check_range("CharPositions::count", l, r, n_);
        const int k = A_.id(c);
        return str_internal::lower_index(pos_.data(), start_[k], start_[k + 1], r) - str_internal::lower_index(pos_.data(), start_[k], start_[k + 1], l);
    }
    template <class C>
    bool has(int l, int r, C c) const { return count(l, r, c) > 0; }
    template <class C>
    int total(C c) const {  // 全体での c の数  O(1)
        const int k = A_.id(c);
        return start_[k + 1] - start_[k];
    }
    template <class C>
    std::vector<int> positions(C c) const {  // c の位置の一覧（昇順）
        const int k = A_.id(c);
        return std::vector<int>(pos_.begin() + start_[k], pos_.begin() + start_[k + 1]);
    }
    int occurrence_index(int i) const {  // s[i] が「その文字の何番目（0 から）の出現」か  O(1)
        if (i < 0 || i >= n_) str_internal::fail("CharPositions::occurrence_index", "位置 " + std::to_string(i) + " は 0〜N-1 の外");
        return occ_[i];
    }

   private:
    void check_pos(const char* func, int i) const {
        if (i < 0 || i > n_) str_internal::fail(func, "位置 " + std::to_string(i) + " は 0〜N（= " + std::to_string(n_) + "）の外");
    }
    Alphabet A_;
    int n_ = 0;
    std::vector<int> start_, pos_, occ_;
};

// ---- ランレングス圧縮：連続する同じ値を {value, len, start} にまとめる ----
template <class T>
struct RunBlock {
    T value{};
    long long len = 0;    // 続く長さ
    long long start = 0;  // 元の列での開始位置
    friend std::ostream& operator<<(std::ostream& os, const RunBlock& b) { return os << "{" << b.value << ", " << b.len << ", " << b.start << "}"; }  // デバッグ用
};
template <class Seq>
std::vector<RunBlock<str_internal::elem_t<Seq>>> run_length(const Seq& s0) {
    const auto& s = str_internal::seq(s0);
    std::vector<RunBlock<str_internal::elem_t<Seq>>> r;
    std::size_t cnt = s.empty() ? 0 : 1;
    for (std::size_t i = 1; i < s.size(); ++i) cnt += !(s[i] == s[i - 1]);
    r.reserve(cnt);
    for (std::size_t i = 0, n = s.size(); i < n;) {
        std::size_t j = i + 1;
        while (j < n && s[j] == s[i]) ++j;
        r.push_back({s[i], (long long)(j - i), (long long)i});
        i = j;
    }
    return r;
}
template <class T>
auto run_length_decode(const std::vector<RunBlock<T>>& runs) {  // 元に戻す（T = char なら string、それ以外は vector<T>）
    using R = std::conditional_t<std::is_same<T, char>::value, std::string, std::vector<T>>;
    long long total = 0;
    for (const auto& b : runs) {
        if (b.len < 0) str_internal::fail("run_length_decode", "len が負");
        total += b.len;
        if (total > 2000000000LL) str_internal::fail("run_length_decode", "元の長さが 2×10^9 を超える");
    }
    R r;
    r.reserve((std::size_t)total);
    for (const auto& b : runs) r.insert(r.end(), (std::size_t)b.len, b.value);
    return r;
}

// ---- 変換 ----
inline std::string to_lower(std::string s) {  // 英大文字を小文字に（ほかはそのまま）
    for (char& c : s)
        if ('A' <= c && c <= 'Z') c = (char)(c - 'A' + 'a');
    return s;
}
inline std::string to_upper(std::string s) {
    for (char& c : s)
        if ('a' <= c && c <= 'z') c = (char)(c - 'a' + 'A');
    return s;
}
inline std::string swap_case(std::string s) {  // 大文字と小文字を入れ替える
    for (char& c : s) {
        if ('A' <= c && c <= 'Z') c = (char)(c - 'A' + 'a');
        else if ('a' <= c && c <= 'z') c = (char)(c - 'a' + 'A');
    }
    return s;
}
template <class Seq>
str_internal::seq_t<Seq> caesar_shift(const Seq& s0, long long k, const Alphabet& A = Alphabet::lower()) {  // A の文字を番号で k ずらす（負も。A にない文字はそのまま）
    auto s = str_internal::seq_t<Seq>(str_internal::seq(s0));
    const long long sg = A.size(), sh = ((k % sg) + sg) % sg;
    for (auto& c : s)
        if (A.contains(c)) c = (std::decay_t<decltype(c)>)A.value((int)((A.id(c) + sh) % sg));
    return s;
}
// 文字の置き換えを重ねる：replace(c, d) で「今 c になっている文字を全部 d に」。apply(s) で s に当てはめる（ABC342 C）
struct CharMap {
    explicit CharMap(const Alphabet& A = Alphabet::lower()) : A_(A), f_(A.size()) {
        for (int i = 0; i < A_.size(); ++i) f_[i] = i;
    }
    template <class C, class D>
    void replace(C c, D d) {  // 今 c になっているものを全部 d に  O(σ)
        const int x = A_.id(c), y = A_.id(d);
        if (x == y) return;
        for (int& v : f_)
            if (v == x) v = y;
    }
    template <class C, class D>
    void set(C c, D d) { f_[A_.id(c)] = A_.id(d); }  // 元の文字 c の行き先を d にする（ほかは変えない）
    template <class C>
    auto operator()(C c) const {  // 元の文字 c の行き先（A にない文字はそのまま）
        return A_.contains(c) ? (C)A_.value(f_[A_.id(c)]) : c;
    }
    template <class Seq>
    str_internal::seq_t<Seq> apply(const Seq& s0) const {  // s の各文字を行き先に（A にない文字はそのまま）
        auto s = str_internal::seq_t<Seq>(str_internal::seq(s0));
        for (auto& c : s) c = (*this)(c);
        return s;
    }
    CharMap compose(const CharMap& g) const {  // 「この置き換えのあとに g」を 1 つにしたもの
        CharMap h(A_);
        for (int i = 0; i < A_.size(); ++i) h.f_[i] = g.f_[f_[i]];
        return h;
    }
    const std::vector<int>& table() const { return f_; }  // 番号 → 行き先の番号

   private:
    Alphabet A_;
    std::vector<int> f_;
};
template <class Seq>
str_internal::seq_t<Seq> reversed(const Seq& s0) {  // 反転したコピー
    const auto& s = str_internal::seq(s0);
    return str_internal::seq_t<Seq>(s.rbegin(), s.rend());
}
template <class Seq>
str_internal::seq_t<Seq> min_with_reverse(const Seq& s0) {  // s と反転のうち辞書順で小さい方（向きを区別しないときの代表。ABC310 C）
    auto s = str_internal::seq_t<Seq>(str_internal::seq(s0));
    auto t = reversed(s);
    return str_internal::lex_compare(t, s) < 0 ? t : s;
}
template <class Seq>
str_internal::seq_t<Seq> sorted_chars(const Seq& s0) {  // 並べ替えたもの（並べ替えて同じかの比較用）
    auto s = str_internal::seq_t<Seq>(str_internal::seq(s0));
    std::sort(s.begin(), s.end(), [](const auto& a, const auto& b) { return str_internal::code_of(a) < str_internal::code_of(b); });
    return s;
}
template <class Seq>
str_internal::seq_t<Seq> repeat_string(const Seq& s0, long long k) {  // s を k 回つなげる
    const auto& s = str_internal::seq(s0);
    if (k < 0) str_internal::fail("repeat_string", "k ≥ 0 で使う");
    if (!s.empty() && k > 2000000000LL / (long long)s.size()) str_internal::fail("repeat_string", "長さが 2×10^9 を超える");
    str_internal::seq_t<Seq> r;
    r.reserve((std::size_t)(k * (long long)s.size()));
    for (long long i = 0; i < k; ++i) r.insert(r.end(), s.begin(), s.end());
    return r;
}

// ---- 比較 ----
template <class A, class B>
int common_prefix_length(const A& a0, const B& b0) {  // 共通の接頭辞の長さ
    const auto& a = str_internal::seq(a0);
    const auto& b = str_internal::seq(b0);
    const std::size_t k = a.size() < b.size() ? a.size() : b.size();
    std::size_t i = 0;
    while (i < k && str_internal::code_of(a[i]) == str_internal::code_of(b[i])) ++i;
    return (int)i;
}
template <class A, class B>
int common_suffix_length(const A& a0, const B& b0) {  // 共通の接尾辞の長さ
    const auto& a = str_internal::seq(a0);
    const auto& b = str_internal::seq(b0);
    const std::size_t k = a.size() < b.size() ? a.size() : b.size();
    std::size_t i = 0;
    while (i < k && str_internal::code_of(a[a.size() - 1 - i]) == str_internal::code_of(b[b.size() - 1 - i])) ++i;
    return (int)i;
}
template <class A, class B>
int hamming_distance(const A& a0, const B& b0) {  // 同じ位置で違う文字の数（長さが同じこと）
    const auto& a = str_internal::seq(a0);
    const auto& b = str_internal::seq(b0);
    if (a.size() != b.size()) str_internal::fail("hamming_distance", "長さが違う");
    int d = 0;
    for (std::size_t i = 0; i < a.size(); ++i) d += str_internal::code_of(a[i]) != str_internal::code_of(b[i]);
    return d;
}
template <class T, class S>
bool is_subsequence(const T& t0, const S& s0) {  // t が s の部分列（飛び飛びでよい）か  O(|s|)。何度も聞くなら 06_subsequence の NextTable
    const auto& t = str_internal::seq(t0);
    const auto& s = str_internal::seq(s0);
    std::size_t j = 0;
    for (std::size_t i = 0; i < s.size() && j < t.size(); ++i)
        if (str_internal::code_of(s[i]) == str_internal::code_of(t[j])) ++j;
    return j == t.size();
}
template <class Seq>
bool is_palindrome(const Seq& s0) {  // 全体が回文か  O(N)。区間を何度も聞くなら 04_palindrome の Palindromes
    const auto& s = str_internal::seq(s0);
    for (std::size_t i = 0, n = s.size(); i < n / 2; ++i)
        if (!(s[i] == s[n - 1 - i])) return false;
    return true;
}

// ---- 区切り・つなげる ----
inline std::vector<std::string> split(const std::string& s, char sep, bool skip_empty = false) {  // sep で区切る（空の部分も残す。skip_empty で除く）
    std::vector<std::string> r;
    std::size_t b = 0;
    for (std::size_t i = 0; i <= s.size(); ++i)
        if (i == s.size() || s[i] == sep) {
            if (!(skip_empty && i == b)) r.push_back(s.substr(b, i - b));
            b = i + 1;
        }
    return r;
}
template <class V>
std::string join(const std::vector<V>& vs, const std::string& sep = " ") {  // 区切りを挟んでつなげる（数も文字列にする。答えの出力に）
    std::string r;
    for (std::size_t i = 0; i < vs.size(); ++i) {
        if (i) r += sep;
        if constexpr (std::is_same<V, std::string>::value) r += vs[i];
        else if constexpr (std::is_same<V, char>::value) r += vs[i];
        else if constexpr (std::is_integral<V>::value) r += std::to_string(vs[i]);
        else {
            std::ostringstream os;
            os << vs[i];
            r += os.str();
        }
    }
    return r;
}
#endif  // GG_STR_BASE
#ifndef GG_SUBSEQUENCE
#define GG_SUBSEQUENCE
// ---- 部分列（文字列第二世代）  次の出現の表・部分列の判定・何文字まで入るか・辞書順最小の部分列・部分列でない最短・部分列の数 ----
// 「部分列」は飛び飛びに選んだもの（順番はそのまま）。連続する「部分文字列」は 02_pattern_match・03_suffix_array
// 数える関数は型 T を受け取る（math_v2 の ModInt998 など。long long なら溢れたら止める）

// ---- 次 / 前の出現を O(1) で（中は 00_str_base の CharPrefix。メモリ (N+1)σ 個の int） ----
struct NextTable {
    NextTable() = default;
    template <class Seq>
    explicit NextTable(const Seq& s, const Alphabet& A = Alphabet::lower()) : cp_(s, A) {}
    int size() const { return cp_.size(); }
    template <class C>
    int next(int i, C c) const { return cp_.next(i, c); }  // i 以上で最初の c（なければ -1。0 ≤ i ≤ N）  O(1)
    template <class C>
    int prev(int i, C c) const { return cp_.prev(i, c); }  // i より前で最後の c（なければ -1）  O(1)
    template <class T>
    int match_length(const T& t0, int from = 0) const {  // s[from..) に t の先頭が何文字まで部分列として入るか（貪欲）  O(|t|)
        const auto& t = str_internal::seq(t0);
        int pos = from, k = 0;
        if (from < 0 || from > size()) str_internal::fail("NextTable::match_length", "from は 0〜N");
        for (; k < (int)t.size(); ++k) {
            const int q = cp_.next(pos, t[k]);
            if (q < 0) break;
            pos = q + 1;
        }
        return k;
    }
    template <class T>
    int match_end(const T& t0, int from = 0) const {  // s[from..) に t を貪欲に当てはめたときの、最後の文字の次の位置（入らなければ -1）  O(|t|)
        const auto& t = str_internal::seq(t0);
        if (from < 0 || from > size()) str_internal::fail("NextTable::match_end", "from は 0〜N");
        int pos = from;
        for (std::size_t k = 0; k < t.size(); ++k) {
            const int q = cp_.next(pos, t[k]);
            if (q < 0) return -1;
            pos = q + 1;
        }
        return pos;
    }
    template <class T>
    bool is_subsequence(const T& t, int from = 0) const { return match_end(t, from) >= 0; }  // t が s[from..) の部分列か  O(|t|)
    const CharPrefix& char_prefix() const { return cp_; }  // 区間の文字の数なども使える

   private:
    CharPrefix cp_;
};

// ---- 先頭から・末尾から何文字まで入るか（ABC324 E） ----
template <class S, class T>
std::vector<int> prefix_match_lengths(const S& s0, const T& t0) {  // r[i] = t の先頭何文字までが s[0, i) の部分列になるか（長さ |s| + 1）  O(|s|)
    const auto& s = str_internal::seq(s0);
    const auto& t = str_internal::seq(t0);
    std::vector<int> r(s.size() + 1, 0);
    int k = 0;
    for (std::size_t i = 0; i < s.size(); ++i) {
        if (k < (int)t.size() && str_internal::code_of(s[i]) == str_internal::code_of(t[k])) ++k;
        r[i + 1] = k;
    }
    return r;
}
template <class S, class T>
std::vector<int> suffix_match_lengths(const S& s0, const T& t0) {  // r[i] = t の末尾何文字までが s[i, N) の部分列になるか（長さ |s| + 1）  O(|s|)
    const auto& s = str_internal::seq(s0);
    const auto& t = str_internal::seq(t0);
    const int n = (int)s.size(), m = (int)t.size();
    std::vector<int> r(n + 1, 0);
    int k = 0;
    for (int i = n - 1; i >= 0; --i) {
        if (k < m && str_internal::code_of(s[i]) == str_internal::code_of(t[m - 1 - k])) ++k;
        r[i] = k;
    }
    return r;
}

// ---- 辞書順最小 / 最大の部分列（長さ k） ----
namespace str_internal {
template <class Seq>
seq_t<Seq> extreme_subsequence(const Seq& s0, int k, bool want_min) {
    const auto& s = seq(s0);
    const int n = (int)s.size();
    if (k < 0 || k > n) fail(want_min ? "lex_smallest_subsequence" : "lex_largest_subsequence", "0 ≤ k ≤ N で使う");
    seq_t<Seq> st;
    for (int i = 0; i < n; ++i) {
        const long long x = code_of(s[i]);
        while (!st.empty() && (want_min ? code_of(st.back()) > x : code_of(st.back()) < x) && (int)st.size() - 1 + (n - i) >= k) st.pop_back();
        if ((int)st.size() < k) st.push_back(s[i]);
    }
    return st;
}
}  // namespace str_internal
template <class Seq>
str_internal::seq_t<Seq> lex_smallest_subsequence(const Seq& s, int k) { return str_internal::extreme_subsequence(s, k, true); }  // 典型90 006  O(N)
template <class Seq>
str_internal::seq_t<Seq> lex_largest_subsequence(const Seq& s, int k) { return str_internal::extreme_subsequence(s, k, false); }

// ---- 部分列にならない最短の文字列（その中で辞書順最小。ARC081 E） ----
template <class Seq>
std::vector<int> shortest_non_subsequence_ids(const Seq& s0, const Alphabet& A = Alphabet::lower()) {  // 番号の列で返す  O(Nσ)
    const auto& s = str_internal::seq(s0);
    const int n = (int)s.size(), sg = A.size();
    std::vector<int> id(n);
    for (int i = 0; i < n; ++i) id[i] = A.id(s[i]);
    // dp[i] = s[i..) の部分列でない最短の長さ。nx[i][c] = i 以上で最初の c（なければ n）
    std::vector<int> dp(n + 2, 0), choice(n + 1, 0), last(sg, n);
    dp[n] = 1;  // 空の残りには長さ 1 の文字列が入らない
    std::vector<int> nxt((std::size_t)(n + 1) * sg);
    for (int c = 0; c < sg; ++c) nxt[(std::size_t)n * sg + c] = n;
    for (int i = n - 1; i >= 0; --i) {
        last[id[i]] = i;
        for (int c = 0; c < sg; ++c) nxt[(std::size_t)i * sg + c] = last[c];
    }
    for (int i = n - 1; i >= 0; --i) {
        int best = 1 << 30, bc = 0;
        for (int c = 0; c < sg; ++c) {
            const int q = nxt[(std::size_t)i * sg + c];
            const int v = q == n ? 1 : 1 + dp[q + 1];
            if (v < best) best = v, bc = c;
        }
        dp[i] = best, choice[i] = bc;
    }
    std::vector<int> r;
    for (int i = 0;;) {
        if (i >= n) {
            r.push_back(0);
            break;
        }
        const int c = choice[i], q = nxt[(std::size_t)i * sg + c];
        r.push_back(c);
        if (q == n) break;
        i = q + 1;
    }
    return r;
}
template <class Seq>
std::string shortest_non_subsequence(const Seq& s, const Alphabet& A = Alphabet::lower()) { return from_ids(shortest_non_subsequence_ids(s, A), A); }  // 文字列で返す

// ---- 部分列の数 ----
template <class T, class S, class P>
T count_subsequence_occurrences(const S& s0, const P& t0) {  // t が s の部分列として現れる数（位置の選び方の数）  O(|s||t|)
    const auto& s = str_internal::seq(s0);
    const auto t = str_internal::codes_of(t0);
    const int m = (int)t.size();
    std::vector<T> dp(m + 1, T(0));  // dp[j] = t の先頭 j 文字の現れ方の数（整数型で溢れたら over[j] に印。答えの dp[m] のときだけ止める）
    std::vector<char> over(m + 1, 0);
    dp[0] = T(1);
    for (std::size_t i = 0; i < s.size(); ++i) {
        const long long x = str_internal::code_of(s[i]);
        for (int j = m; j >= 1; --j)
            if (t[j - 1] == x) {
                bool o = over[j];
                str_internal::add_count(dp[j], o, dp[j - 1], over[j - 1]);
                over[j] = o;
            }
    }
    if (over[m]) str_internal::count_overflow("count_subsequence_occurrences");
    return dp[m];
}
template <class T, class Seq>
T count_distinct_subsequences(const Seq& s0, bool include_empty = false, int min_step = 1) {  // 異なる部分列の数（選ぶ位置の間を min_step 以上あける。ABC214 F は 2）  O(N log N)（値の圧縮）
    if (min_step < 1) str_internal::fail("count_distinct_subsequences", "min_step ≥ 1 で使う");
    const CompressedValues cv = compress_values(s0);
    const int n = (int)cv.ids.size(), g = min_step;
    std::vector<int> last(cv.size(), -1);
    std::vector<T> f(n, T(0)), P(n + 1, T(0));  // f[i] = 貪欲に左から当てはめると i で終わる異なる部分列の数、P = f の累積和
    T total = include_empty ? T(1) : T(0);
    for (int i = 0; i < n; ++i) {
        const int prev = last[cv.ids[i]];
        last[cv.ids[i]] = i;
        T v = prev < 0 ? T(1) : T(0);
        const int lo = prev < 0 ? 0 : std::max(0, prev - g + 1), hi = i - g;  // 1 つ前に選んだ位置 j の範囲 [lo, hi]
        if (lo <= hi) v = str_internal::checked_add(v, P[hi + 1] - P[lo], "count_distinct_subsequences");
        f[i] = v;
        P[i + 1] = str_internal::checked_add(P[i], v, "count_distinct_subsequences");
        total = str_internal::checked_add(total, v, "count_distinct_subsequences");
    }
    return total;
}
template <class T, class Seq>
std::vector<T> distinct_subsequences_by_length(const Seq& s0) {  // r[len] = 長さ len の異なる部分列の数（len = 0〜N。r[0] = 1）  O(N^2)
    const CompressedValues cv = compress_values(s0);
    const int n = (int)cv.ids.size();
    std::vector<int> prev(n), last(cv.size(), -1);
    for (int i = 0; i < n; ++i) prev[i] = last[cv.ids[i]], last[cv.ids[i]] = i;
    std::vector<T> r(n + 1, T(0)), f(n, T(0)), P(n + 1, T(0));
    r[0] = T(1);
    for (int i = 0; i < n; ++i) f[i] = prev[i] < 0 ? T(1) : T(0);  // 長さ 1
    for (int len = 1; len <= n; ++len) {
        T sum = T(0);
        for (int i = 0; i < n; ++i) sum = str_internal::checked_add(sum, f[i], "distinct_subsequences_by_length");
        r[len] = sum;
        if (len == n) break;
        P[0] = T(0);
        for (int i = 0; i < n; ++i) P[i + 1] = str_internal::checked_add(P[i], f[i], "distinct_subsequences_by_length");
        for (int i = 0; i < n; ++i) {  // 次の長さ：j ∈ [max(0, prev), i - 1]
            const int lo = std::max(0, prev[i]), hi = i - 1;
            f[i] = lo <= hi ? P[hi + 1] - P[lo] : T(0);
        }
    }
    return r;
}
#endif  // GG_SUBSEQUENCE

#ifndef GG_STR_BASE
#define GG_STR_BASE
// ---- 文字列基盤（文字列第二世代）  文字の種類 Alphabet・文字の数と位置・ランレングス・変換・比較・区切り ----
// 列（string / vector<整数>）を受け取る関数はどちらでもよい。文字列リテラル（"abc"）も渡せる
// 範囲は半開区間 [l, r)（0 始まり）。位置が「ない」ときは -1。使える範囲の外の引数はメッセージを出して止める
namespace str_internal {
[[noreturn]] inline void fail(const char* func, const std::string& msg) {
    std::cerr << func << ": " << msg << "\n";
    std::abort();
}
// 文字列リテラルを std::string に（それ以外はそのまま参照を返す）
template <class T>
const T& seq(const T& s) { return s; }
template <std::size_t N>
std::string seq(const char (&s)[N]) { return std::string(s); }
inline std::string seq(const char* s) { return std::string(s); }
template <class T>
using seq_t = std::decay_t<decltype(seq(std::declval<const T&>()))>;
template <class T>
using elem_t = std::decay_t<decltype(std::declval<const seq_t<T>&>()[0])>;
// 要素を「値」に（char は unsigned char として 0〜255。負の添字にしない）
template <class T>
long long code_of(T x) {
    if constexpr (std::is_same<T, char>::value || std::is_same<T, signed char>::value || std::is_same<T, unsigned char>::value)
        return (long long)(unsigned char)x;
    else
        return (long long)x;
}
// 範囲の確認（止めるときの文字列作りは別の関数に分けて、確認そのものを軽くする）
[[noreturn]] __attribute__((noinline, cold)) inline void range_fail(const char* func, long long l, long long r, long long n) {
    fail(func, "範囲 [l, r) は 0 ≤ l ≤ r ≤ 長さ で指定する（l = " + std::to_string(l) + ", r = " + std::to_string(r) + ", 長さ = " + std::to_string(n) + "）");
}
inline void check_range(const char* func, long long l, long long r, long long n) {
    if (__builtin_expect(l < 0 || r > n || l > r, 0)) range_fail(func, l, r, n);
}
// a[lo, hi) の中で x 以上の最初の位置（a は昇順）。分岐の少ない二分探索
inline int lower_index(const int* a, int lo, int hi, int x) {
    int len = hi - lo;
    if (len <= 0) return lo;
    const int* b = a + lo;
    while (len > 1) {
        const int half = len >> 1;
        b += (b[half - 1] < x) ? half : 0;
        len -= half;
    }
    return (int)(b - a) + (*b < x);
}
// 昇順の a の中で x より大きい最初の位置（std::upper_bound と同じ。デバッグビルドでも O(log N)）
template <class T>
int upper_index(const std::vector<T>& a, const T& x) {
    int lo = 0, hi = (int)a.size();
    while (lo < hi) {
        const int mid = (int)((unsigned)(lo + hi) >> 1);
        if (!(x < a[mid])) lo = mid + 1;
        else hi = mid;
    }
    return lo;
}
// 数え上げの足し算：T が整数型なら溢れたら止める（mint・double などはそのまま）
template <class T>
T checked_add(const T& a, const T& b, const char* func) {
    if constexpr (std::is_integral<T>::value) {
        T r;
        if (__builtin_add_overflow(a, b, &r)) fail(func, "値が型に収まらない（答えを mod で求めるなら math_v2 の ModInt998 などを使う）");
        return r;
    } else {
        return a + b;
    }
}
// 数え上げの途中の足し算：T が整数型なら溢れたことを印 oa に残すだけ（途中で溢れても答えに使わなければ止めない。最後に印を見て止める）
template <class T>
void add_count(T& a, bool& oa, const T& b, bool ob) {
    if constexpr (std::is_integral<T>::value) {
        T r;
        oa = __builtin_add_overflow(a, b, &r) || oa || ob;
        a = r;
    } else {
        a = a + b;
    }
}
[[noreturn]] inline void count_overflow(const char* func) { fail(func, "値が型に収まらない（答えを mod で求めるなら math_v2 の ModInt998 などを使う）"); }
// 列 → 値の列（char は 0〜255）
template <class Seq>
std::vector<long long> codes_of(const Seq& s0) {
    const auto& s = seq(s0);
    std::vector<long long> c(s.size());
    for (std::size_t i = 0; i < s.size(); ++i) c[i] = code_of(s[i]);
    return c;
}
// 辞書順の比較（-1 / 0 / 1）。char は unsigned char として比べる（std::string と同じ順）
template <class A, class B>
int lex_compare(const A& a, const B& b) {
    const std::size_t n = a.size(), m = b.size(), k = n < m ? n : m;
    for (std::size_t i = 0; i < k; ++i) {
        const long long x = code_of(a[i]), y = code_of(b[i]);
        if (x != y) return x < y ? -1 : 1;
    }
    return n == m ? 0 : (n < m ? -1 : 1);
}
}  // namespace str_internal

// ---- 部分文字列の位置 [l, r) ----
struct StrSpan {
    int l = 0, r = 0;
    int len() const { return r - l; }
    bool empty() const { return r <= l; }
    friend bool operator==(const StrSpan& a, const StrSpan& b) { return a.l == b.l && a.r == b.r; }
    friend bool operator!=(const StrSpan& a, const StrSpan& b) { return !(a == b); }
    friend std::ostream& operator<<(std::ostream& os, const StrSpan& sp) { return os << "[" << sp.l << ", " << sp.r << ")"; }  // デバッグ用
};
template <class Seq>
str_internal::seq_t<Seq> substr_of(const Seq& s0, StrSpan sp) {  // s[sp.l, sp.r) のコピー（string なら string、vector なら vector）
    const auto& s = str_internal::seq(s0);
    str_internal::check_range("substr_of", sp.l, sp.r, (long long)s.size());
    return str_internal::seq_t<Seq>(s.begin() + sp.l, s.begin() + sp.r);
}

// ---- 文字の種類（Alphabet）：文字・値 → 番号 0〜size()-1。番号の順が辞書順として扱われる ----
//   Alphabet::lower()（a〜z。既定）/ upper()（A〜Z）/ digits()（0〜9）/ letters()（A〜Z, a〜z の 52 種。ASCII 順）
//   alnum()（0〜9, A〜Z, a〜z の 62 種。ASCII 順）/ ascii()（0〜127）/ bytes()（0〜255。日本語など 128 以上の char も）
//   of("()[]")（並べた文字だけ。並べた順が番号の順 = 辞書順）/ ints(K)（整数 0〜K-1）/ ints(lo, hi)（整数 [lo, hi)）
//   範囲外の文字・値を id に渡すとメッセージを出して止まる（contains で先に確かめられる）
struct Alphabet {
    Alphabet() : Alphabet(of("abcdefghijklmnopqrstuvwxyz")) {}
    static Alphabet lower() { return of("abcdefghijklmnopqrstuvwxyz"); }
    static Alphabet upper() { return of("ABCDEFGHIJKLMNOPQRSTUVWXYZ"); }
    static Alphabet digits() { return of("0123456789"); }
    static Alphabet letters() { return of("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"); }
    static Alphabet alnum() { return of("0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"); }
    static Alphabet ascii() { return code_range(0, 128); }
    static Alphabet bytes() { return code_range(0, 256); }
    static Alphabet of(const std::string& cs) {  // 並べた文字（同じ文字を 2 回並べると止まる）
        if (cs.empty()) str_internal::fail("Alphabet::of", "文字が 1 つもない");
        Alphabet A(0);
        A.kind_ = CHARS;
        for (std::size_t i = 0; i < cs.size(); ++i) {
            const unsigned char c = (unsigned char)cs[i];
            if (A.table_[c] >= 0) str_internal::fail("Alphabet::of", std::string("同じ文字が 2 回ある: '") + cs[i] + "'");
            A.table_[c] = (short)i;
            A.chars_[i] = c;
        }
        A.size_ = (int)cs.size();
        return A;
    }
    static Alphabet ints(long long K) { return ints(0, K); }  // 整数 0〜K-1
    static Alphabet ints(long long lo, long long hi) {         // 整数 [lo, hi)
        if (hi <= lo || (unsigned long long)hi - (unsigned long long)lo > 1000000000ULL)
            str_internal::fail("Alphabet::ints", "lo < hi かつ hi - lo ≤ 10^9 で指定する（大きな値は compress_values で 0〜K-1 に）");
        Alphabet A(0);
        A.kind_ = INTS;
        A.lo_ = lo;
        A.size_ = (int)(hi - lo);
        return A;
    }

    int size() const { return size_; }
    bool is_chars() const { return kind_ == CHARS; }  // 文字の表で決めたものか（of・lower など）。false なら ints
    template <class T>
    bool contains(T x) const { return index_of(str_internal::code_of(x)) >= 0; }
    template <class T>
    int id(T x) const {  // 文字・値の番号（範囲外は止める）
        const long long v = str_internal::code_of(x);
        const int r = index_of(v);
        if (r < 0) out_of_range(v);
        return r;
    }
    long long value(int i) const {  // 番号 i の値（文字なら文字コード 0〜255）
        check_index(i);
        return kind_ == CHARS ? (long long)chars_[i] : lo_ + i;
    }
    char ch(int i) const {  // 番号 i の文字（ints で 0〜255 の外の値なら止める）
        const long long v = value(i);
        if (v < 0 || v > 255) str_internal::fail("Alphabet::ch", "番号 " + std::to_string(i) + " の値 " + std::to_string(v) + " は文字にできない");
        return (char)(unsigned char)v;
    }
    std::string chars() const {  // 文字を番号の順に並べたもの（ints なら 0〜255 の範囲のときだけ）
        std::string s;
        for (int i = 0; i < size_; ++i) s += ch(i);
        return s;
    }

   private:
    enum Kind : unsigned char { CHARS, INTS };
    explicit Alphabet(int) { table_.fill(-1); }
    static Alphabet code_range(int lo, int hi) {
        Alphabet A(0);
        A.kind_ = CHARS;
        for (int c = lo; c < hi; ++c) A.table_[c] = (short)(c - lo), A.chars_[c - lo] = (unsigned char)c;
        A.size_ = hi - lo;
        return A;
    }
    int index_of(long long v) const {
        if (kind_ == CHARS) return (unsigned long long)v < 256ULL ? table_[(std::size_t)v] : -1;
        const unsigned long long d = (unsigned long long)v - (unsigned long long)lo_;
        return d < (unsigned long long)size_ ? (int)d : -1;
    }
    void check_index(int i) const {
        if (i < 0 || i >= size_) str_internal::fail("Alphabet", "番号 " + std::to_string(i) + " は 0〜" + std::to_string(size_ - 1) + " の外");
    }
    [[noreturn]] void out_of_range(long long v) const {
        std::string what = "値 " + std::to_string(v) + " は";
        if (kind_ == CHARS && v >= 32 && v < 127) what = std::string("文字 '") + (char)v + "'（コード " + std::to_string(v) + "）は";
        str_internal::fail("Alphabet::id", what + "文字の種類に入っていない（種類の数 " + std::to_string(size_) + "）");
    }
    Kind kind_ = CHARS;
    int size_ = 0;
    long long lo_ = 0;
    std::array<short, 256> table_{};
    std::array<unsigned char, 256> chars_{};
};

template <class Seq>
std::vector<int> to_ids(const Seq& s0, const Alphabet& A = Alphabet::lower()) {  // 各文字の番号の列
    const auto& s = str_internal::seq(s0);
    std::vector<int> r(s.size());
    for (std::size_t i = 0; i < s.size(); ++i) r[i] = A.id(s[i]);
    return r;
}
inline std::string from_ids(const std::vector<int>& ids, const Alphabet& A = Alphabet::lower()) {  // 番号の列 → 文字列
    std::string s(ids.size(), '\0');
    for (std::size_t i = 0; i < ids.size(); ++i) s[i] = A.ch(ids[i]);
    return s;
}

// ---- 数の列を 0〜K-1 に圧縮（大小の順を保つ）：Alphabet::ints(K) と組み合わせる ----
struct CompressedValues {
    std::vector<int> ids;            // 各要素の番号（0〜K-1。値の小さい順）
    std::vector<long long> values;   // 番号 → 元の値（昇順。K = values.size()）
    int size() const { return (int)values.size(); }
    Alphabet alphabet() const { return Alphabet::ints(values.empty() ? 1 : (long long)values.size()); }
};
template <class Seq>
CompressedValues compress_values(const Seq& s0) {
    const auto& s = str_internal::seq(s0);
    const std::size_t n = s.size();
    std::vector<std::pair<long long, int>> vi(n);
    for (std::size_t i = 0; i < n; ++i) vi[i] = {str_internal::code_of(s[i]), (int)i};
    std::sort(vi.begin(), vi.end());
    CompressedValues c;
    c.ids.resize(n);
    for (std::size_t i = 0; i < n; ++i) {
        if (i == 0 || vi[i].first != vi[i - 1].first) c.values.push_back(vi[i].first);
        c.ids[vi[i].second] = (int)c.values.size() - 1;
    }
    return c;
}

// ---- 文字の数 ----
template <class Seq>
std::vector<int> count_chars(const Seq& s0, const Alphabet& A = Alphabet::lower()) {  // 番号ごとの個数（長さ A.size()）
    const auto& s = str_internal::seq(s0);
    std::vector<int> c(A.size(), 0);
    for (std::size_t i = 0; i < s.size(); ++i) ++c[A.id(s[i])];
    return c;
}

// ---- 区間の文字の数・次 / 前の出現を O(1) で（前処理 O(Nσ)・メモリ (N+1)σ 個の int。N = 2×10^5, σ = 26 で約 21MB） ----
//   メモリを抑えたいとき・σ が大きいときは CharPositions（O(log N)・メモリ O(N + σ)）
struct CharPrefix {
    CharPrefix() = default;
    template <class Seq>
    explicit CharPrefix(const Seq& s0, const Alphabet& A = Alphabet::lower()) : A_(A) {
        const auto& s = str_internal::seq(s0);
        n_ = (int)s.size(), sg_ = A_.size();
        if ((long long)(n_ + 1) * sg_ > 600000000LL) str_internal::fail("CharPrefix", "(N + 1)σ が 6×10^8 を超える（メモリが足りない。CharPositions を使う）");
        cnt_.assign((std::size_t)(n_ + 1) * sg_, 0);
        std::vector<int> id(n_);
        start_.assign(sg_ + 1, 0);
        for (int i = 0; i < n_; ++i) {
            const int* src = cnt_.data() + (std::size_t)i * sg_;
            int* dst = cnt_.data() + (std::size_t)(i + 1) * sg_;
            std::copy(src, src + sg_, dst);
            id[i] = A_.id(s[i]);
            ++dst[id[i]], ++start_[id[i] + 1];
        }
        for (int k = 0; k < sg_; ++k) start_[k + 1] += start_[k];
        pos_.resize(n_);
        std::vector<int> cur(start_.begin(), start_.end() - 1);
        for (int i = 0; i < n_; ++i) pos_[cur[id[i]]++] = i;
    }
    int size() const { return n_; }
    const Alphabet& alphabet() const { return A_; }
    template <class C>
    int count(int l, int r, C c) const {  // [l, r) にある c の数  O(1)
        str_internal::check_range("CharPrefix::count", l, r, n_);
        const int k = A_.id(c);
        return cnt_[(std::size_t)r * sg_ + k] - cnt_[(std::size_t)l * sg_ + k];
    }
    template <class C>
    bool has(int l, int r, C c) const { return count(l, r, c) > 0; }  // [l, r) に c があるか  O(1)
    std::vector<int> counts(int l, int r) const {  // [l, r) の番号ごとの数  O(σ)
        str_internal::check_range("CharPrefix::counts", l, r, n_);
        std::vector<int> c(sg_);
        for (int k = 0; k < sg_; ++k) c[k] = cnt_[(std::size_t)r * sg_ + k] - cnt_[(std::size_t)l * sg_ + k];
        return c;
    }
    int distinct(int l, int r) const {  // [l, r) にある文字の種類の数  O(σ)
        str_internal::check_range("CharPrefix::distinct", l, r, n_);
        int d = 0;
        for (int k = 0; k < sg_; ++k) d += cnt_[(std::size_t)r * sg_ + k] != cnt_[(std::size_t)l * sg_ + k];
        return d;
    }
    template <class C>
    int next(int i, C c) const {  // i 以上で最初の c の位置（なければ -1。0 ≤ i ≤ N）  O(1)
        check_pos("CharPrefix::next", i);
        const int k = A_.id(c), m = cnt_[(std::size_t)i * sg_ + k];
        return start_[k] + m < start_[k + 1] ? pos_[start_[k] + m] : -1;
    }
    template <class C>
    int prev(int i, C c) const {  // i より前（[0, i)）で最後の c の位置（なければ -1。0 ≤ i ≤ N）  O(1)
        check_pos("CharPrefix::prev", i);
        const int k = A_.id(c), m = cnt_[(std::size_t)i * sg_ + k];
        return m > 0 ? pos_[start_[k] + m - 1] : -1;
    }
    template <class C>
    int kth(C c, long long k) const {  // c の k 番目（0 から）の位置（なければ -1）  O(1)
        const int id = A_.id(c);
        return (k >= 0 && k < start_[id + 1] - start_[id]) ? pos_[start_[id] + k] : -1;
    }
    template <class C>
    int total(C c) const {  // 全体での c の数  O(1)
        const int k = A_.id(c);
        return start_[k + 1] - start_[k];
    }

   private:
    void check_pos(const char* func, int i) const {
        if (i < 0 || i > n_) str_internal::fail(func, "位置 " + std::to_string(i) + " は 0〜N（= " + std::to_string(n_) + "）の外");
    }
    Alphabet A_;
    int n_ = 0, sg_ = 0;
    std::vector<int> cnt_;           // cnt_[i * σ + k] = [0, i) にある番号 k の数
    std::vector<int> start_, pos_;   // 番号ごとの出現位置（pos_[start_[k] .. start_[k+1])）
};

// ---- 文字の出現位置（前処理 O(N + σ)・メモリ O(N + σ)。σ が大きくても使える） ----
struct CharPositions {
    CharPositions() = default;
    template <class Seq>
    explicit CharPositions(const Seq& s0, const Alphabet& A = Alphabet::lower()) : A_(A) {
        const auto& s = str_internal::seq(s0);
        n_ = (int)s.size();
        const int sg = A_.size();
        std::vector<int> id(n_);
        start_.assign(sg + 1, 0);
        for (int i = 0; i < n_; ++i) id[i] = A_.id(s[i]), ++start_[id[i] + 1];
        for (int k = 0; k < sg; ++k) start_[k + 1] += start_[k];
        pos_.resize(n_), occ_.resize(n_);
        std::vector<int> cur(start_.begin(), start_.end() - 1);
        for (int i = 0; i < n_; ++i) occ_[i] = cur[id[i]] - start_[id[i]], pos_[cur[id[i]]++] = i;
    }
    int size() const { return n_; }
    const Alphabet& alphabet() const { return A_; }
    template <class C>
    int next(int i, C c) const {  // i 以上で最初の c の位置（なければ -1。0 ≤ i ≤ N）  O(log N)
        check_pos("CharPositions::next", i);
        const int k = A_.id(c), j = str_internal::lower_index(pos_.data(), start_[k], start_[k + 1], i);
        return j < start_[k + 1] ? pos_[j] : -1;
    }
    template <class C>
    int prev(int i, C c) const {  // i より前（[0, i)）で最後の c の位置（なければ -1。0 ≤ i ≤ N）  O(log N)
        check_pos("CharPositions::prev", i);
        const int k = A_.id(c), j = str_internal::lower_index(pos_.data(), start_[k], start_[k + 1], i);
        return j > start_[k] ? pos_[j - 1] : -1;
    }
    template <class C>
    int kth(C c, long long k) const {  // c の k 番目（0 から）の位置（なければ -1）  O(1)
        const int id = A_.id(c);
        return (k >= 0 && k < start_[id + 1] - start_[id]) ? pos_[start_[id] + k] : -1;
    }
    template <class C>
    int count(int l, int r, C c) const {  // [l, r) にある c の数  O(log N)
        str_internal::check_range("CharPositions::count", l, r, n_);
        const int k = A_.id(c);
        return str_internal::lower_index(pos_.data(), start_[k], start_[k + 1], r) - str_internal::lower_index(pos_.data(), start_[k], start_[k + 1], l);
    }
    template <class C>
    bool has(int l, int r, C c) const { return count(l, r, c) > 0; }
    template <class C>
    int total(C c) const {  // 全体での c の数  O(1)
        const int k = A_.id(c);
        return start_[k + 1] - start_[k];
    }
    template <class C>
    std::vector<int> positions(C c) const {  // c の位置の一覧（昇順）
        const int k = A_.id(c);
        return std::vector<int>(pos_.begin() + start_[k], pos_.begin() + start_[k + 1]);
    }
    int occurrence_index(int i) const {  // s[i] が「その文字の何番目（0 から）の出現」か  O(1)
        if (i < 0 || i >= n_) str_internal::fail("CharPositions::occurrence_index", "位置 " + std::to_string(i) + " は 0〜N-1 の外");
        return occ_[i];
    }

   private:
    void check_pos(const char* func, int i) const {
        if (i < 0 || i > n_) str_internal::fail(func, "位置 " + std::to_string(i) + " は 0〜N（= " + std::to_string(n_) + "）の外");
    }
    Alphabet A_;
    int n_ = 0;
    std::vector<int> start_, pos_, occ_;
};

// ---- ランレングス圧縮：連続する同じ値を {value, len, start} にまとめる ----
template <class T>
struct RunBlock {
    T value{};
    long long len = 0;    // 続く長さ
    long long start = 0;  // 元の列での開始位置
    friend std::ostream& operator<<(std::ostream& os, const RunBlock& b) { return os << "{" << b.value << ", " << b.len << ", " << b.start << "}"; }  // デバッグ用
};
template <class Seq>
std::vector<RunBlock<str_internal::elem_t<Seq>>> run_length(const Seq& s0) {
    const auto& s = str_internal::seq(s0);
    std::vector<RunBlock<str_internal::elem_t<Seq>>> r;
    std::size_t cnt = s.empty() ? 0 : 1;
    for (std::size_t i = 1; i < s.size(); ++i) cnt += !(s[i] == s[i - 1]);
    r.reserve(cnt);
    for (std::size_t i = 0, n = s.size(); i < n;) {
        std::size_t j = i + 1;
        while (j < n && s[j] == s[i]) ++j;
        r.push_back({s[i], (long long)(j - i), (long long)i});
        i = j;
    }
    return r;
}
template <class T>
auto run_length_decode(const std::vector<RunBlock<T>>& runs) {  // 元に戻す（T = char なら string、それ以外は vector<T>）
    using R = std::conditional_t<std::is_same<T, char>::value, std::string, std::vector<T>>;
    long long total = 0;
    for (const auto& b : runs) {
        if (b.len < 0) str_internal::fail("run_length_decode", "len が負");
        total += b.len;
        if (total > 2000000000LL) str_internal::fail("run_length_decode", "元の長さが 2×10^9 を超える");
    }
    R r;
    r.reserve((std::size_t)total);
    for (const auto& b : runs) r.insert(r.end(), (std::size_t)b.len, b.value);
    return r;
}

// ---- 変換 ----
inline std::string to_lower(std::string s) {  // 英大文字を小文字に（ほかはそのまま）
    for (char& c : s)
        if ('A' <= c && c <= 'Z') c = (char)(c - 'A' + 'a');
    return s;
}
inline std::string to_upper(std::string s) {
    for (char& c : s)
        if ('a' <= c && c <= 'z') c = (char)(c - 'a' + 'A');
    return s;
}
inline std::string swap_case(std::string s) {  // 大文字と小文字を入れ替える
    for (char& c : s) {
        if ('A' <= c && c <= 'Z') c = (char)(c - 'A' + 'a');
        else if ('a' <= c && c <= 'z') c = (char)(c - 'a' + 'A');
    }
    return s;
}
template <class Seq>
str_internal::seq_t<Seq> caesar_shift(const Seq& s0, long long k, const Alphabet& A = Alphabet::lower()) {  // A の文字を番号で k ずらす（負も。A にない文字はそのまま）
    auto s = str_internal::seq_t<Seq>(str_internal::seq(s0));
    const long long sg = A.size(), sh = ((k % sg) + sg) % sg;
    for (auto& c : s)
        if (A.contains(c)) c = (std::decay_t<decltype(c)>)A.value((int)((A.id(c) + sh) % sg));
    return s;
}
// 文字の置き換えを重ねる：replace(c, d) で「今 c になっている文字を全部 d に」。apply(s) で s に当てはめる（ABC342 C）
struct CharMap {
    explicit CharMap(const Alphabet& A = Alphabet::lower()) : A_(A), f_(A.size()) {
        for (int i = 0; i < A_.size(); ++i) f_[i] = i;
    }
    template <class C, class D>
    void replace(C c, D d) {  // 今 c になっているものを全部 d に  O(σ)
        const int x = A_.id(c), y = A_.id(d);
        if (x == y) return;
        for (int& v : f_)
            if (v == x) v = y;
    }
    template <class C, class D>
    void set(C c, D d) { f_[A_.id(c)] = A_.id(d); }  // 元の文字 c の行き先を d にする（ほかは変えない）
    template <class C>
    auto operator()(C c) const {  // 元の文字 c の行き先（A にない文字はそのまま）
        return A_.contains(c) ? (C)A_.value(f_[A_.id(c)]) : c;
    }
    template <class Seq>
    str_internal::seq_t<Seq> apply(const Seq& s0) const {  // s の各文字を行き先に（A にない文字はそのまま）
        auto s = str_internal::seq_t<Seq>(str_internal::seq(s0));
        for (auto& c : s) c = (*this)(c);
        return s;
    }
    CharMap compose(const CharMap& g) const {  // 「この置き換えのあとに g」を 1 つにしたもの
        CharMap h(A_);
        for (int i = 0; i < A_.size(); ++i) h.f_[i] = g.f_[f_[i]];
        return h;
    }
    const std::vector<int>& table() const { return f_; }  // 番号 → 行き先の番号

   private:
    Alphabet A_;
    std::vector<int> f_;
};
template <class Seq>
str_internal::seq_t<Seq> reversed(const Seq& s0) {  // 反転したコピー
    const auto& s = str_internal::seq(s0);
    return str_internal::seq_t<Seq>(s.rbegin(), s.rend());
}
template <class Seq>
str_internal::seq_t<Seq> min_with_reverse(const Seq& s0) {  // s と反転のうち辞書順で小さい方（向きを区別しないときの代表。ABC310 C）
    auto s = str_internal::seq_t<Seq>(str_internal::seq(s0));
    auto t = reversed(s);
    return str_internal::lex_compare(t, s) < 0 ? t : s;
}
template <class Seq>
str_internal::seq_t<Seq> sorted_chars(const Seq& s0) {  // 並べ替えたもの（並べ替えて同じかの比較用）
    auto s = str_internal::seq_t<Seq>(str_internal::seq(s0));
    std::sort(s.begin(), s.end(), [](const auto& a, const auto& b) { return str_internal::code_of(a) < str_internal::code_of(b); });
    return s;
}
template <class Seq>
str_internal::seq_t<Seq> repeat_string(const Seq& s0, long long k) {  // s を k 回つなげる
    const auto& s = str_internal::seq(s0);
    if (k < 0) str_internal::fail("repeat_string", "k ≥ 0 で使う");
    if (!s.empty() && k > 2000000000LL / (long long)s.size()) str_internal::fail("repeat_string", "長さが 2×10^9 を超える");
    str_internal::seq_t<Seq> r;
    r.reserve((std::size_t)(k * (long long)s.size()));
    for (long long i = 0; i < k; ++i) r.insert(r.end(), s.begin(), s.end());
    return r;
}

// ---- 比較 ----
template <class A, class B>
int common_prefix_length(const A& a0, const B& b0) {  // 共通の接頭辞の長さ
    const auto& a = str_internal::seq(a0);
    const auto& b = str_internal::seq(b0);
    const std::size_t k = a.size() < b.size() ? a.size() : b.size();
    std::size_t i = 0;
    while (i < k && str_internal::code_of(a[i]) == str_internal::code_of(b[i])) ++i;
    return (int)i;
}
template <class A, class B>
int common_suffix_length(const A& a0, const B& b0) {  // 共通の接尾辞の長さ
    const auto& a = str_internal::seq(a0);
    const auto& b = str_internal::seq(b0);
    const std::size_t k = a.size() < b.size() ? a.size() : b.size();
    std::size_t i = 0;
    while (i < k && str_internal::code_of(a[a.size() - 1 - i]) == str_internal::code_of(b[b.size() - 1 - i])) ++i;
    return (int)i;
}
template <class A, class B>
int hamming_distance(const A& a0, const B& b0) {  // 同じ位置で違う文字の数（長さが同じこと）
    const auto& a = str_internal::seq(a0);
    const auto& b = str_internal::seq(b0);
    if (a.size() != b.size()) str_internal::fail("hamming_distance", "長さが違う");
    int d = 0;
    for (std::size_t i = 0; i < a.size(); ++i) d += str_internal::code_of(a[i]) != str_internal::code_of(b[i]);
    return d;
}
template <class T, class S>
bool is_subsequence(const T& t0, const S& s0) {  // t が s の部分列（飛び飛びでよい）か  O(|s|)。何度も聞くなら 06_subsequence の NextTable
    const auto& t = str_internal::seq(t0);
    const auto& s = str_internal::seq(s0);
    std::size_t j = 0;
    for (std::size_t i = 0; i < s.size() && j < t.size(); ++i)
        if (str_internal::code_of(s[i]) == str_internal::code_of(t[j])) ++j;
    return j == t.size();
}
template <class Seq>
bool is_palindrome(const Seq& s0) {  // 全体が回文か  O(N)。区間を何度も聞くなら 04_palindrome の Palindromes
    const auto& s = str_internal::seq(s0);
    for (std::size_t i = 0, n = s.size(); i < n / 2; ++i)
        if (!(s[i] == s[n - 1 - i])) return false;
    return true;
}

// ---- 区切り・つなげる ----
inline std::vector<std::string> split(const std::string& s, char sep, bool skip_empty = false) {  // sep で区切る（空の部分も残す。skip_empty で除く）
    std::vector<std::string> r;
    std::size_t b = 0;
    for (std::size_t i = 0; i <= s.size(); ++i)
        if (i == s.size() || s[i] == sep) {
            if (!(skip_empty && i == b)) r.push_back(s.substr(b, i - b));
            b = i + 1;
        }
    return r;
}
template <class V>
std::string join(const std::vector<V>& vs, const std::string& sep = " ") {  // 区切りを挟んでつなげる（数も文字列にする。答えの出力に）
    std::string r;
    for (std::size_t i = 0; i < vs.size(); ++i) {
        if (i) r += sep;
        if constexpr (std::is_same<V, std::string>::value) r += vs[i];
        else if constexpr (std::is_same<V, char>::value) r += vs[i];
        else if constexpr (std::is_integral<V>::value) r += std::to_string(vs[i]);
        else {
            std::ostringstream os;
            os << vs[i];
            r += os.str();
        }
    }
    return r;
}
#endif  // GG_STR_BASE
#ifndef GG_BRACKET
#define GG_BRACKET
// ---- 括弧列（文字列第二世代） ----
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

#ifndef GG_STR_BASE
#define GG_STR_BASE
// ---- 文字列基盤（文字列第二世代）  文字の種類 Alphabet・文字の数と位置・ランレングス・変換・比較・区切り ----
// 列（string / vector<整数>）を受け取る関数はどちらでもよい。文字列リテラル（"abc"）も渡せる
// 範囲は半開区間 [l, r)（0 始まり）。位置が「ない」ときは -1。使える範囲の外の引数はメッセージを出して止める
namespace str_internal {
[[noreturn]] inline void fail(const char* func, const std::string& msg) {
    std::cerr << func << ": " << msg << "\n";
    std::abort();
}
// 文字列リテラルを std::string に（それ以外はそのまま参照を返す）
template <class T>
const T& seq(const T& s) { return s; }
template <std::size_t N>
std::string seq(const char (&s)[N]) { return std::string(s); }
inline std::string seq(const char* s) { return std::string(s); }
template <class T>
using seq_t = std::decay_t<decltype(seq(std::declval<const T&>()))>;
template <class T>
using elem_t = std::decay_t<decltype(std::declval<const seq_t<T>&>()[0])>;
// 要素を「値」に（char は unsigned char として 0〜255。負の添字にしない）
template <class T>
long long code_of(T x) {
    if constexpr (std::is_same<T, char>::value || std::is_same<T, signed char>::value || std::is_same<T, unsigned char>::value)
        return (long long)(unsigned char)x;
    else
        return (long long)x;
}
// 範囲の確認（止めるときの文字列作りは別の関数に分けて、確認そのものを軽くする）
[[noreturn]] __attribute__((noinline, cold)) inline void range_fail(const char* func, long long l, long long r, long long n) {
    fail(func, "範囲 [l, r) は 0 ≤ l ≤ r ≤ 長さ で指定する（l = " + std::to_string(l) + ", r = " + std::to_string(r) + ", 長さ = " + std::to_string(n) + "）");
}
inline void check_range(const char* func, long long l, long long r, long long n) {
    if (__builtin_expect(l < 0 || r > n || l > r, 0)) range_fail(func, l, r, n);
}
// a[lo, hi) の中で x 以上の最初の位置（a は昇順）。分岐の少ない二分探索
inline int lower_index(const int* a, int lo, int hi, int x) {
    int len = hi - lo;
    if (len <= 0) return lo;
    const int* b = a + lo;
    while (len > 1) {
        const int half = len >> 1;
        b += (b[half - 1] < x) ? half : 0;
        len -= half;
    }
    return (int)(b - a) + (*b < x);
}
// 昇順の a の中で x より大きい最初の位置（std::upper_bound と同じ。デバッグビルドでも O(log N)）
template <class T>
int upper_index(const std::vector<T>& a, const T& x) {
    int lo = 0, hi = (int)a.size();
    while (lo < hi) {
        const int mid = (int)((unsigned)(lo + hi) >> 1);
        if (!(x < a[mid])) lo = mid + 1;
        else hi = mid;
    }
    return lo;
}
// 数え上げの足し算：T が整数型なら溢れたら止める（mint・double などはそのまま）
template <class T>
T checked_add(const T& a, const T& b, const char* func) {
    if constexpr (std::is_integral<T>::value) {
        T r;
        if (__builtin_add_overflow(a, b, &r)) fail(func, "値が型に収まらない（答えを mod で求めるなら math_v2 の ModInt998 などを使う）");
        return r;
    } else {
        return a + b;
    }
}
// 数え上げの途中の足し算：T が整数型なら溢れたことを印 oa に残すだけ（途中で溢れても答えに使わなければ止めない。最後に印を見て止める）
template <class T>
void add_count(T& a, bool& oa, const T& b, bool ob) {
    if constexpr (std::is_integral<T>::value) {
        T r;
        oa = __builtin_add_overflow(a, b, &r) || oa || ob;
        a = r;
    } else {
        a = a + b;
    }
}
[[noreturn]] inline void count_overflow(const char* func) { fail(func, "値が型に収まらない（答えを mod で求めるなら math_v2 の ModInt998 などを使う）"); }
// 列 → 値の列（char は 0〜255）
template <class Seq>
std::vector<long long> codes_of(const Seq& s0) {
    const auto& s = seq(s0);
    std::vector<long long> c(s.size());
    for (std::size_t i = 0; i < s.size(); ++i) c[i] = code_of(s[i]);
    return c;
}
// 辞書順の比較（-1 / 0 / 1）。char は unsigned char として比べる（std::string と同じ順）
template <class A, class B>
int lex_compare(const A& a, const B& b) {
    const std::size_t n = a.size(), m = b.size(), k = n < m ? n : m;
    for (std::size_t i = 0; i < k; ++i) {
        const long long x = code_of(a[i]), y = code_of(b[i]);
        if (x != y) return x < y ? -1 : 1;
    }
    return n == m ? 0 : (n < m ? -1 : 1);
}
}  // namespace str_internal

// ---- 部分文字列の位置 [l, r) ----
struct StrSpan {
    int l = 0, r = 0;
    int len() const { return r - l; }
    bool empty() const { return r <= l; }
    friend bool operator==(const StrSpan& a, const StrSpan& b) { return a.l == b.l && a.r == b.r; }
    friend bool operator!=(const StrSpan& a, const StrSpan& b) { return !(a == b); }
    friend std::ostream& operator<<(std::ostream& os, const StrSpan& sp) { return os << "[" << sp.l << ", " << sp.r << ")"; }  // デバッグ用
};
template <class Seq>
str_internal::seq_t<Seq> substr_of(const Seq& s0, StrSpan sp) {  // s[sp.l, sp.r) のコピー（string なら string、vector なら vector）
    const auto& s = str_internal::seq(s0);
    str_internal::check_range("substr_of", sp.l, sp.r, (long long)s.size());
    return str_internal::seq_t<Seq>(s.begin() + sp.l, s.begin() + sp.r);
}

// ---- 文字の種類（Alphabet）：文字・値 → 番号 0〜size()-1。番号の順が辞書順として扱われる ----
//   Alphabet::lower()（a〜z。既定）/ upper()（A〜Z）/ digits()（0〜9）/ letters()（A〜Z, a〜z の 52 種。ASCII 順）
//   alnum()（0〜9, A〜Z, a〜z の 62 種。ASCII 順）/ ascii()（0〜127）/ bytes()（0〜255。日本語など 128 以上の char も）
//   of("()[]")（並べた文字だけ。並べた順が番号の順 = 辞書順）/ ints(K)（整数 0〜K-1）/ ints(lo, hi)（整数 [lo, hi)）
//   範囲外の文字・値を id に渡すとメッセージを出して止まる（contains で先に確かめられる）
struct Alphabet {
    Alphabet() : Alphabet(of("abcdefghijklmnopqrstuvwxyz")) {}
    static Alphabet lower() { return of("abcdefghijklmnopqrstuvwxyz"); }
    static Alphabet upper() { return of("ABCDEFGHIJKLMNOPQRSTUVWXYZ"); }
    static Alphabet digits() { return of("0123456789"); }
    static Alphabet letters() { return of("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"); }
    static Alphabet alnum() { return of("0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"); }
    static Alphabet ascii() { return code_range(0, 128); }
    static Alphabet bytes() { return code_range(0, 256); }
    static Alphabet of(const std::string& cs) {  // 並べた文字（同じ文字を 2 回並べると止まる）
        if (cs.empty()) str_internal::fail("Alphabet::of", "文字が 1 つもない");
        Alphabet A(0);
        A.kind_ = CHARS;
        for (std::size_t i = 0; i < cs.size(); ++i) {
            const unsigned char c = (unsigned char)cs[i];
            if (A.table_[c] >= 0) str_internal::fail("Alphabet::of", std::string("同じ文字が 2 回ある: '") + cs[i] + "'");
            A.table_[c] = (short)i;
            A.chars_[i] = c;
        }
        A.size_ = (int)cs.size();
        return A;
    }
    static Alphabet ints(long long K) { return ints(0, K); }  // 整数 0〜K-1
    static Alphabet ints(long long lo, long long hi) {         // 整数 [lo, hi)
        if (hi <= lo || (unsigned long long)hi - (unsigned long long)lo > 1000000000ULL)
            str_internal::fail("Alphabet::ints", "lo < hi かつ hi - lo ≤ 10^9 で指定する（大きな値は compress_values で 0〜K-1 に）");
        Alphabet A(0);
        A.kind_ = INTS;
        A.lo_ = lo;
        A.size_ = (int)(hi - lo);
        return A;
    }

    int size() const { return size_; }
    bool is_chars() const { return kind_ == CHARS; }  // 文字の表で決めたものか（of・lower など）。false なら ints
    template <class T>
    bool contains(T x) const { return index_of(str_internal::code_of(x)) >= 0; }
    template <class T>
    int id(T x) const {  // 文字・値の番号（範囲外は止める）
        const long long v = str_internal::code_of(x);
        const int r = index_of(v);
        if (r < 0) out_of_range(v);
        return r;
    }
    long long value(int i) const {  // 番号 i の値（文字なら文字コード 0〜255）
        check_index(i);
        return kind_ == CHARS ? (long long)chars_[i] : lo_ + i;
    }
    char ch(int i) const {  // 番号 i の文字（ints で 0〜255 の外の値なら止める）
        const long long v = value(i);
        if (v < 0 || v > 255) str_internal::fail("Alphabet::ch", "番号 " + std::to_string(i) + " の値 " + std::to_string(v) + " は文字にできない");
        return (char)(unsigned char)v;
    }
    std::string chars() const {  // 文字を番号の順に並べたもの（ints なら 0〜255 の範囲のときだけ）
        std::string s;
        for (int i = 0; i < size_; ++i) s += ch(i);
        return s;
    }

   private:
    enum Kind : unsigned char { CHARS, INTS };
    explicit Alphabet(int) { table_.fill(-1); }
    static Alphabet code_range(int lo, int hi) {
        Alphabet A(0);
        A.kind_ = CHARS;
        for (int c = lo; c < hi; ++c) A.table_[c] = (short)(c - lo), A.chars_[c - lo] = (unsigned char)c;
        A.size_ = hi - lo;
        return A;
    }
    int index_of(long long v) const {
        if (kind_ == CHARS) return (unsigned long long)v < 256ULL ? table_[(std::size_t)v] : -1;
        const unsigned long long d = (unsigned long long)v - (unsigned long long)lo_;
        return d < (unsigned long long)size_ ? (int)d : -1;
    }
    void check_index(int i) const {
        if (i < 0 || i >= size_) str_internal::fail("Alphabet", "番号 " + std::to_string(i) + " は 0〜" + std::to_string(size_ - 1) + " の外");
    }
    [[noreturn]] void out_of_range(long long v) const {
        std::string what = "値 " + std::to_string(v) + " は";
        if (kind_ == CHARS && v >= 32 && v < 127) what = std::string("文字 '") + (char)v + "'（コード " + std::to_string(v) + "）は";
        str_internal::fail("Alphabet::id", what + "文字の種類に入っていない（種類の数 " + std::to_string(size_) + "）");
    }
    Kind kind_ = CHARS;
    int size_ = 0;
    long long lo_ = 0;
    std::array<short, 256> table_{};
    std::array<unsigned char, 256> chars_{};
};

template <class Seq>
std::vector<int> to_ids(const Seq& s0, const Alphabet& A = Alphabet::lower()) {  // 各文字の番号の列
    const auto& s = str_internal::seq(s0);
    std::vector<int> r(s.size());
    for (std::size_t i = 0; i < s.size(); ++i) r[i] = A.id(s[i]);
    return r;
}
inline std::string from_ids(const std::vector<int>& ids, const Alphabet& A = Alphabet::lower()) {  // 番号の列 → 文字列
    std::string s(ids.size(), '\0');
    for (std::size_t i = 0; i < ids.size(); ++i) s[i] = A.ch(ids[i]);
    return s;
}

// ---- 数の列を 0〜K-1 に圧縮（大小の順を保つ）：Alphabet::ints(K) と組み合わせる ----
struct CompressedValues {
    std::vector<int> ids;            // 各要素の番号（0〜K-1。値の小さい順）
    std::vector<long long> values;   // 番号 → 元の値（昇順。K = values.size()）
    int size() const { return (int)values.size(); }
    Alphabet alphabet() const { return Alphabet::ints(values.empty() ? 1 : (long long)values.size()); }
};
template <class Seq>
CompressedValues compress_values(const Seq& s0) {
    const auto& s = str_internal::seq(s0);
    const std::size_t n = s.size();
    std::vector<std::pair<long long, int>> vi(n);
    for (std::size_t i = 0; i < n; ++i) vi[i] = {str_internal::code_of(s[i]), (int)i};
    std::sort(vi.begin(), vi.end());
    CompressedValues c;
    c.ids.resize(n);
    for (std::size_t i = 0; i < n; ++i) {
        if (i == 0 || vi[i].first != vi[i - 1].first) c.values.push_back(vi[i].first);
        c.ids[vi[i].second] = (int)c.values.size() - 1;
    }
    return c;
}

// ---- 文字の数 ----
template <class Seq>
std::vector<int> count_chars(const Seq& s0, const Alphabet& A = Alphabet::lower()) {  // 番号ごとの個数（長さ A.size()）
    const auto& s = str_internal::seq(s0);
    std::vector<int> c(A.size(), 0);
    for (std::size_t i = 0; i < s.size(); ++i) ++c[A.id(s[i])];
    return c;
}

// ---- 区間の文字の数・次 / 前の出現を O(1) で（前処理 O(Nσ)・メモリ (N+1)σ 個の int。N = 2×10^5, σ = 26 で約 21MB） ----
//   メモリを抑えたいとき・σ が大きいときは CharPositions（O(log N)・メモリ O(N + σ)）
struct CharPrefix {
    CharPrefix() = default;
    template <class Seq>
    explicit CharPrefix(const Seq& s0, const Alphabet& A = Alphabet::lower()) : A_(A) {
        const auto& s = str_internal::seq(s0);
        n_ = (int)s.size(), sg_ = A_.size();
        if ((long long)(n_ + 1) * sg_ > 600000000LL) str_internal::fail("CharPrefix", "(N + 1)σ が 6×10^8 を超える（メモリが足りない。CharPositions を使う）");
        cnt_.assign((std::size_t)(n_ + 1) * sg_, 0);
        std::vector<int> id(n_);
        start_.assign(sg_ + 1, 0);
        for (int i = 0; i < n_; ++i) {
            const int* src = cnt_.data() + (std::size_t)i * sg_;
            int* dst = cnt_.data() + (std::size_t)(i + 1) * sg_;
            std::copy(src, src + sg_, dst);
            id[i] = A_.id(s[i]);
            ++dst[id[i]], ++start_[id[i] + 1];
        }
        for (int k = 0; k < sg_; ++k) start_[k + 1] += start_[k];
        pos_.resize(n_);
        std::vector<int> cur(start_.begin(), start_.end() - 1);
        for (int i = 0; i < n_; ++i) pos_[cur[id[i]]++] = i;
    }
    int size() const { return n_; }
    const Alphabet& alphabet() const { return A_; }
    template <class C>
    int count(int l, int r, C c) const {  // [l, r) にある c の数  O(1)
        str_internal::check_range("CharPrefix::count", l, r, n_);
        const int k = A_.id(c);
        return cnt_[(std::size_t)r * sg_ + k] - cnt_[(std::size_t)l * sg_ + k];
    }
    template <class C>
    bool has(int l, int r, C c) const { return count(l, r, c) > 0; }  // [l, r) に c があるか  O(1)
    std::vector<int> counts(int l, int r) const {  // [l, r) の番号ごとの数  O(σ)
        str_internal::check_range("CharPrefix::counts", l, r, n_);
        std::vector<int> c(sg_);
        for (int k = 0; k < sg_; ++k) c[k] = cnt_[(std::size_t)r * sg_ + k] - cnt_[(std::size_t)l * sg_ + k];
        return c;
    }
    int distinct(int l, int r) const {  // [l, r) にある文字の種類の数  O(σ)
        str_internal::check_range("CharPrefix::distinct", l, r, n_);
        int d = 0;
        for (int k = 0; k < sg_; ++k) d += cnt_[(std::size_t)r * sg_ + k] != cnt_[(std::size_t)l * sg_ + k];
        return d;
    }
    template <class C>
    int next(int i, C c) const {  // i 以上で最初の c の位置（なければ -1。0 ≤ i ≤ N）  O(1)
        check_pos("CharPrefix::next", i);
        const int k = A_.id(c), m = cnt_[(std::size_t)i * sg_ + k];
        return start_[k] + m < start_[k + 1] ? pos_[start_[k] + m] : -1;
    }
    template <class C>
    int prev(int i, C c) const {  // i より前（[0, i)）で最後の c の位置（なければ -1。0 ≤ i ≤ N）  O(1)
        check_pos("CharPrefix::prev", i);
        const int k = A_.id(c), m = cnt_[(std::size_t)i * sg_ + k];
        return m > 0 ? pos_[start_[k] + m - 1] : -1;
    }
    template <class C>
    int kth(C c, long long k) const {  // c の k 番目（0 から）の位置（なければ -1）  O(1)
        const int id = A_.id(c);
        return (k >= 0 && k < start_[id + 1] - start_[id]) ? pos_[start_[id] + k] : -1;
    }
    template <class C>
    int total(C c) const {  // 全体での c の数  O(1)
        const int k = A_.id(c);
        return start_[k + 1] - start_[k];
    }

   private:
    void check_pos(const char* func, int i) const {
        if (i < 0 || i > n_) str_internal::fail(func, "位置 " + std::to_string(i) + " は 0〜N（= " + std::to_string(n_) + "）の外");
    }
    Alphabet A_;
    int n_ = 0, sg_ = 0;
    std::vector<int> cnt_;           // cnt_[i * σ + k] = [0, i) にある番号 k の数
    std::vector<int> start_, pos_;   // 番号ごとの出現位置（pos_[start_[k] .. start_[k+1])）
};

// ---- 文字の出現位置（前処理 O(N + σ)・メモリ O(N + σ)。σ が大きくても使える） ----
struct CharPositions {
    CharPositions() = default;
    template <class Seq>
    explicit CharPositions(const Seq& s0, const Alphabet& A = Alphabet::lower()) : A_(A) {
        const auto& s = str_internal::seq(s0);
        n_ = (int)s.size();
        const int sg = A_.size();
        std::vector<int> id(n_);
        start_.assign(sg + 1, 0);
        for (int i = 0; i < n_; ++i) id[i] = A_.id(s[i]), ++start_[id[i] + 1];
        for (int k = 0; k < sg; ++k) start_[k + 1] += start_[k];
        pos_.resize(n_), occ_.resize(n_);
        std::vector<int> cur(start_.begin(), start_.end() - 1);
        for (int i = 0; i < n_; ++i) occ_[i] = cur[id[i]] - start_[id[i]], pos_[cur[id[i]]++] = i;
    }
    int size() const { return n_; }
    const Alphabet& alphabet() const { return A_; }
    template <class C>
    int next(int i, C c) const {  // i 以上で最初の c の位置（なければ -1。0 ≤ i ≤ N）  O(log N)
        check_pos("CharPositions::next", i);
        const int k = A_.id(c), j = str_internal::lower_index(pos_.data(), start_[k], start_[k + 1], i);
        return j < start_[k + 1] ? pos_[j] : -1;
    }
    template <class C>
    int prev(int i, C c) const {  // i より前（[0, i)）で最後の c の位置（なければ -1。0 ≤ i ≤ N）  O(log N)
        check_pos("CharPositions::prev", i);
        const int k = A_.id(c), j = str_internal::lower_index(pos_.data(), start_[k], start_[k + 1], i);
        return j > start_[k] ? pos_[j - 1] : -1;
    }
    template <class C>
    int kth(C c, long long k) const {  // c の k 番目（0 から）の位置（なければ -1）  O(1)
        const int id = A_.id(c);
        return (k >= 0 && k < start_[id + 1] - start_[id]) ? pos_[start_[id] + k] : -1;
    }
    template <class C>
    int count(int l, int r, C c) const {  // [l, r) にある c の数  O(log N)
        str_internal::check_range("CharPositions::count", l, r, n_);
        const int k = A_.id(c);
        return str_internal::lower_index(pos_.data(), start_[k], start_[k + 1], r) - str_internal::lower_index(pos_.data(), start_[k], start_[k + 1], l);
    }
    template <class C>
    bool has(int l, int r, C c) const { return count(l, r, c) > 0; }
    template <class C>
    int total(C c) const {  // 全体での c の数  O(1)
        const int k = A_.id(c);
        return start_[k + 1] - start_[k];
    }
    template <class C>
    std::vector<int> positions(C c) const {  // c の位置の一覧（昇順）
        const int k = A_.id(c);
        return std::vector<int>(pos_.begin() + start_[k], pos_.begin() + start_[k + 1]);
    }
    int occurrence_index(int i) const {  // s[i] が「その文字の何番目（0 から）の出現」か  O(1)
        if (i < 0 || i >= n_) str_internal::fail("CharPositions::occurrence_index", "位置 " + std::to_string(i) + " は 0〜N-1 の外");
        return occ_[i];
    }

   private:
    void check_pos(const char* func, int i) const {
        if (i < 0 || i > n_) str_internal::fail(func, "位置 " + std::to_string(i) + " は 0〜N（= " + std::to_string(n_) + "）の外");
    }
    Alphabet A_;
    int n_ = 0;
    std::vector<int> start_, pos_, occ_;
};

// ---- ランレングス圧縮：連続する同じ値を {value, len, start} にまとめる ----
template <class T>
struct RunBlock {
    T value{};
    long long len = 0;    // 続く長さ
    long long start = 0;  // 元の列での開始位置
    friend std::ostream& operator<<(std::ostream& os, const RunBlock& b) { return os << "{" << b.value << ", " << b.len << ", " << b.start << "}"; }  // デバッグ用
};
template <class Seq>
std::vector<RunBlock<str_internal::elem_t<Seq>>> run_length(const Seq& s0) {
    const auto& s = str_internal::seq(s0);
    std::vector<RunBlock<str_internal::elem_t<Seq>>> r;
    std::size_t cnt = s.empty() ? 0 : 1;
    for (std::size_t i = 1; i < s.size(); ++i) cnt += !(s[i] == s[i - 1]);
    r.reserve(cnt);
    for (std::size_t i = 0, n = s.size(); i < n;) {
        std::size_t j = i + 1;
        while (j < n && s[j] == s[i]) ++j;
        r.push_back({s[i], (long long)(j - i), (long long)i});
        i = j;
    }
    return r;
}
template <class T>
auto run_length_decode(const std::vector<RunBlock<T>>& runs) {  // 元に戻す（T = char なら string、それ以外は vector<T>）
    using R = std::conditional_t<std::is_same<T, char>::value, std::string, std::vector<T>>;
    long long total = 0;
    for (const auto& b : runs) {
        if (b.len < 0) str_internal::fail("run_length_decode", "len が負");
        total += b.len;
        if (total > 2000000000LL) str_internal::fail("run_length_decode", "元の長さが 2×10^9 を超える");
    }
    R r;
    r.reserve((std::size_t)total);
    for (const auto& b : runs) r.insert(r.end(), (std::size_t)b.len, b.value);
    return r;
}

// ---- 変換 ----
inline std::string to_lower(std::string s) {  // 英大文字を小文字に（ほかはそのまま）
    for (char& c : s)
        if ('A' <= c && c <= 'Z') c = (char)(c - 'A' + 'a');
    return s;
}
inline std::string to_upper(std::string s) {
    for (char& c : s)
        if ('a' <= c && c <= 'z') c = (char)(c - 'a' + 'A');
    return s;
}
inline std::string swap_case(std::string s) {  // 大文字と小文字を入れ替える
    for (char& c : s) {
        if ('A' <= c && c <= 'Z') c = (char)(c - 'A' + 'a');
        else if ('a' <= c && c <= 'z') c = (char)(c - 'a' + 'A');
    }
    return s;
}
template <class Seq>
str_internal::seq_t<Seq> caesar_shift(const Seq& s0, long long k, const Alphabet& A = Alphabet::lower()) {  // A の文字を番号で k ずらす（負も。A にない文字はそのまま）
    auto s = str_internal::seq_t<Seq>(str_internal::seq(s0));
    const long long sg = A.size(), sh = ((k % sg) + sg) % sg;
    for (auto& c : s)
        if (A.contains(c)) c = (std::decay_t<decltype(c)>)A.value((int)((A.id(c) + sh) % sg));
    return s;
}
// 文字の置き換えを重ねる：replace(c, d) で「今 c になっている文字を全部 d に」。apply(s) で s に当てはめる（ABC342 C）
struct CharMap {
    explicit CharMap(const Alphabet& A = Alphabet::lower()) : A_(A), f_(A.size()) {
        for (int i = 0; i < A_.size(); ++i) f_[i] = i;
    }
    template <class C, class D>
    void replace(C c, D d) {  // 今 c になっているものを全部 d に  O(σ)
        const int x = A_.id(c), y = A_.id(d);
        if (x == y) return;
        for (int& v : f_)
            if (v == x) v = y;
    }
    template <class C, class D>
    void set(C c, D d) { f_[A_.id(c)] = A_.id(d); }  // 元の文字 c の行き先を d にする（ほかは変えない）
    template <class C>
    auto operator()(C c) const {  // 元の文字 c の行き先（A にない文字はそのまま）
        return A_.contains(c) ? (C)A_.value(f_[A_.id(c)]) : c;
    }
    template <class Seq>
    str_internal::seq_t<Seq> apply(const Seq& s0) const {  // s の各文字を行き先に（A にない文字はそのまま）
        auto s = str_internal::seq_t<Seq>(str_internal::seq(s0));
        for (auto& c : s) c = (*this)(c);
        return s;
    }
    CharMap compose(const CharMap& g) const {  // 「この置き換えのあとに g」を 1 つにしたもの
        CharMap h(A_);
        for (int i = 0; i < A_.size(); ++i) h.f_[i] = g.f_[f_[i]];
        return h;
    }
    const std::vector<int>& table() const { return f_; }  // 番号 → 行き先の番号

   private:
    Alphabet A_;
    std::vector<int> f_;
};
template <class Seq>
str_internal::seq_t<Seq> reversed(const Seq& s0) {  // 反転したコピー
    const auto& s = str_internal::seq(s0);
    return str_internal::seq_t<Seq>(s.rbegin(), s.rend());
}
template <class Seq>
str_internal::seq_t<Seq> min_with_reverse(const Seq& s0) {  // s と反転のうち辞書順で小さい方（向きを区別しないときの代表。ABC310 C）
    auto s = str_internal::seq_t<Seq>(str_internal::seq(s0));
    auto t = reversed(s);
    return str_internal::lex_compare(t, s) < 0 ? t : s;
}
template <class Seq>
str_internal::seq_t<Seq> sorted_chars(const Seq& s0) {  // 並べ替えたもの（並べ替えて同じかの比較用）
    auto s = str_internal::seq_t<Seq>(str_internal::seq(s0));
    std::sort(s.begin(), s.end(), [](const auto& a, const auto& b) { return str_internal::code_of(a) < str_internal::code_of(b); });
    return s;
}
template <class Seq>
str_internal::seq_t<Seq> repeat_string(const Seq& s0, long long k) {  // s を k 回つなげる
    const auto& s = str_internal::seq(s0);
    if (k < 0) str_internal::fail("repeat_string", "k ≥ 0 で使う");
    if (!s.empty() && k > 2000000000LL / (long long)s.size()) str_internal::fail("repeat_string", "長さが 2×10^9 を超える");
    str_internal::seq_t<Seq> r;
    r.reserve((std::size_t)(k * (long long)s.size()));
    for (long long i = 0; i < k; ++i) r.insert(r.end(), s.begin(), s.end());
    return r;
}

// ---- 比較 ----
template <class A, class B>
int common_prefix_length(const A& a0, const B& b0) {  // 共通の接頭辞の長さ
    const auto& a = str_internal::seq(a0);
    const auto& b = str_internal::seq(b0);
    const std::size_t k = a.size() < b.size() ? a.size() : b.size();
    std::size_t i = 0;
    while (i < k && str_internal::code_of(a[i]) == str_internal::code_of(b[i])) ++i;
    return (int)i;
}
template <class A, class B>
int common_suffix_length(const A& a0, const B& b0) {  // 共通の接尾辞の長さ
    const auto& a = str_internal::seq(a0);
    const auto& b = str_internal::seq(b0);
    const std::size_t k = a.size() < b.size() ? a.size() : b.size();
    std::size_t i = 0;
    while (i < k && str_internal::code_of(a[a.size() - 1 - i]) == str_internal::code_of(b[b.size() - 1 - i])) ++i;
    return (int)i;
}
template <class A, class B>
int hamming_distance(const A& a0, const B& b0) {  // 同じ位置で違う文字の数（長さが同じこと）
    const auto& a = str_internal::seq(a0);
    const auto& b = str_internal::seq(b0);
    if (a.size() != b.size()) str_internal::fail("hamming_distance", "長さが違う");
    int d = 0;
    for (std::size_t i = 0; i < a.size(); ++i) d += str_internal::code_of(a[i]) != str_internal::code_of(b[i]);
    return d;
}
template <class T, class S>
bool is_subsequence(const T& t0, const S& s0) {  // t が s の部分列（飛び飛びでよい）か  O(|s|)。何度も聞くなら 06_subsequence の NextTable
    const auto& t = str_internal::seq(t0);
    const auto& s = str_internal::seq(s0);
    std::size_t j = 0;
    for (std::size_t i = 0; i < s.size() && j < t.size(); ++i)
        if (str_internal::code_of(s[i]) == str_internal::code_of(t[j])) ++j;
    return j == t.size();
}
template <class Seq>
bool is_palindrome(const Seq& s0) {  // 全体が回文か  O(N)。区間を何度も聞くなら 04_palindrome の Palindromes
    const auto& s = str_internal::seq(s0);
    for (std::size_t i = 0, n = s.size(); i < n / 2; ++i)
        if (!(s[i] == s[n - 1 - i])) return false;
    return true;
}

// ---- 区切り・つなげる ----
inline std::vector<std::string> split(const std::string& s, char sep, bool skip_empty = false) {  // sep で区切る（空の部分も残す。skip_empty で除く）
    std::vector<std::string> r;
    std::size_t b = 0;
    for (std::size_t i = 0; i <= s.size(); ++i)
        if (i == s.size() || s[i] == sep) {
            if (!(skip_empty && i == b)) r.push_back(s.substr(b, i - b));
            b = i + 1;
        }
    return r;
}
template <class V>
std::string join(const std::vector<V>& vs, const std::string& sep = " ") {  // 区切りを挟んでつなげる（数も文字列にする。答えの出力に）
    std::string r;
    for (std::size_t i = 0; i < vs.size(); ++i) {
        if (i) r += sep;
        if constexpr (std::is_same<V, std::string>::value) r += vs[i];
        else if constexpr (std::is_same<V, char>::value) r += vs[i];
        else if constexpr (std::is_integral<V>::value) r += std::to_string(vs[i]);
        else {
            std::ostringstream os;
            os << vs[i];
            r += os.str();
        }
    }
    return r;
}
#endif  // GG_STR_BASE
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

#ifndef GG_STR_BASE
#define GG_STR_BASE
// ---- 文字列基盤（文字列第二世代）  文字の種類 Alphabet・文字の数と位置・ランレングス・変換・比較・区切り ----
// 列（string / vector<整数>）を受け取る関数はどちらでもよい。文字列リテラル（"abc"）も渡せる
// 範囲は半開区間 [l, r)（0 始まり）。位置が「ない」ときは -1。使える範囲の外の引数はメッセージを出して止める
namespace str_internal {
[[noreturn]] inline void fail(const char* func, const std::string& msg) {
    std::cerr << func << ": " << msg << "\n";
    std::abort();
}
// 文字列リテラルを std::string に（それ以外はそのまま参照を返す）
template <class T>
const T& seq(const T& s) { return s; }
template <std::size_t N>
std::string seq(const char (&s)[N]) { return std::string(s); }
inline std::string seq(const char* s) { return std::string(s); }
template <class T>
using seq_t = std::decay_t<decltype(seq(std::declval<const T&>()))>;
template <class T>
using elem_t = std::decay_t<decltype(std::declval<const seq_t<T>&>()[0])>;
// 要素を「値」に（char は unsigned char として 0〜255。負の添字にしない）
template <class T>
long long code_of(T x) {
    if constexpr (std::is_same<T, char>::value || std::is_same<T, signed char>::value || std::is_same<T, unsigned char>::value)
        return (long long)(unsigned char)x;
    else
        return (long long)x;
}
// 範囲の確認（止めるときの文字列作りは別の関数に分けて、確認そのものを軽くする）
[[noreturn]] __attribute__((noinline, cold)) inline void range_fail(const char* func, long long l, long long r, long long n) {
    fail(func, "範囲 [l, r) は 0 ≤ l ≤ r ≤ 長さ で指定する（l = " + std::to_string(l) + ", r = " + std::to_string(r) + ", 長さ = " + std::to_string(n) + "）");
}
inline void check_range(const char* func, long long l, long long r, long long n) {
    if (__builtin_expect(l < 0 || r > n || l > r, 0)) range_fail(func, l, r, n);
}
// a[lo, hi) の中で x 以上の最初の位置（a は昇順）。分岐の少ない二分探索
inline int lower_index(const int* a, int lo, int hi, int x) {
    int len = hi - lo;
    if (len <= 0) return lo;
    const int* b = a + lo;
    while (len > 1) {
        const int half = len >> 1;
        b += (b[half - 1] < x) ? half : 0;
        len -= half;
    }
    return (int)(b - a) + (*b < x);
}
// 昇順の a の中で x より大きい最初の位置（std::upper_bound と同じ。デバッグビルドでも O(log N)）
template <class T>
int upper_index(const std::vector<T>& a, const T& x) {
    int lo = 0, hi = (int)a.size();
    while (lo < hi) {
        const int mid = (int)((unsigned)(lo + hi) >> 1);
        if (!(x < a[mid])) lo = mid + 1;
        else hi = mid;
    }
    return lo;
}
// 数え上げの足し算：T が整数型なら溢れたら止める（mint・double などはそのまま）
template <class T>
T checked_add(const T& a, const T& b, const char* func) {
    if constexpr (std::is_integral<T>::value) {
        T r;
        if (__builtin_add_overflow(a, b, &r)) fail(func, "値が型に収まらない（答えを mod で求めるなら math_v2 の ModInt998 などを使う）");
        return r;
    } else {
        return a + b;
    }
}
// 数え上げの途中の足し算：T が整数型なら溢れたことを印 oa に残すだけ（途中で溢れても答えに使わなければ止めない。最後に印を見て止める）
template <class T>
void add_count(T& a, bool& oa, const T& b, bool ob) {
    if constexpr (std::is_integral<T>::value) {
        T r;
        oa = __builtin_add_overflow(a, b, &r) || oa || ob;
        a = r;
    } else {
        a = a + b;
    }
}
[[noreturn]] inline void count_overflow(const char* func) { fail(func, "値が型に収まらない（答えを mod で求めるなら math_v2 の ModInt998 などを使う）"); }
// 列 → 値の列（char は 0〜255）
template <class Seq>
std::vector<long long> codes_of(const Seq& s0) {
    const auto& s = seq(s0);
    std::vector<long long> c(s.size());
    for (std::size_t i = 0; i < s.size(); ++i) c[i] = code_of(s[i]);
    return c;
}
// 辞書順の比較（-1 / 0 / 1）。char は unsigned char として比べる（std::string と同じ順）
template <class A, class B>
int lex_compare(const A& a, const B& b) {
    const std::size_t n = a.size(), m = b.size(), k = n < m ? n : m;
    for (std::size_t i = 0; i < k; ++i) {
        const long long x = code_of(a[i]), y = code_of(b[i]);
        if (x != y) return x < y ? -1 : 1;
    }
    return n == m ? 0 : (n < m ? -1 : 1);
}
}  // namespace str_internal

// ---- 部分文字列の位置 [l, r) ----
struct StrSpan {
    int l = 0, r = 0;
    int len() const { return r - l; }
    bool empty() const { return r <= l; }
    friend bool operator==(const StrSpan& a, const StrSpan& b) { return a.l == b.l && a.r == b.r; }
    friend bool operator!=(const StrSpan& a, const StrSpan& b) { return !(a == b); }
    friend std::ostream& operator<<(std::ostream& os, const StrSpan& sp) { return os << "[" << sp.l << ", " << sp.r << ")"; }  // デバッグ用
};
template <class Seq>
str_internal::seq_t<Seq> substr_of(const Seq& s0, StrSpan sp) {  // s[sp.l, sp.r) のコピー（string なら string、vector なら vector）
    const auto& s = str_internal::seq(s0);
    str_internal::check_range("substr_of", sp.l, sp.r, (long long)s.size());
    return str_internal::seq_t<Seq>(s.begin() + sp.l, s.begin() + sp.r);
}

// ---- 文字の種類（Alphabet）：文字・値 → 番号 0〜size()-1。番号の順が辞書順として扱われる ----
//   Alphabet::lower()（a〜z。既定）/ upper()（A〜Z）/ digits()（0〜9）/ letters()（A〜Z, a〜z の 52 種。ASCII 順）
//   alnum()（0〜9, A〜Z, a〜z の 62 種。ASCII 順）/ ascii()（0〜127）/ bytes()（0〜255。日本語など 128 以上の char も）
//   of("()[]")（並べた文字だけ。並べた順が番号の順 = 辞書順）/ ints(K)（整数 0〜K-1）/ ints(lo, hi)（整数 [lo, hi)）
//   範囲外の文字・値を id に渡すとメッセージを出して止まる（contains で先に確かめられる）
struct Alphabet {
    Alphabet() : Alphabet(of("abcdefghijklmnopqrstuvwxyz")) {}
    static Alphabet lower() { return of("abcdefghijklmnopqrstuvwxyz"); }
    static Alphabet upper() { return of("ABCDEFGHIJKLMNOPQRSTUVWXYZ"); }
    static Alphabet digits() { return of("0123456789"); }
    static Alphabet letters() { return of("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"); }
    static Alphabet alnum() { return of("0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"); }
    static Alphabet ascii() { return code_range(0, 128); }
    static Alphabet bytes() { return code_range(0, 256); }
    static Alphabet of(const std::string& cs) {  // 並べた文字（同じ文字を 2 回並べると止まる）
        if (cs.empty()) str_internal::fail("Alphabet::of", "文字が 1 つもない");
        Alphabet A(0);
        A.kind_ = CHARS;
        for (std::size_t i = 0; i < cs.size(); ++i) {
            const unsigned char c = (unsigned char)cs[i];
            if (A.table_[c] >= 0) str_internal::fail("Alphabet::of", std::string("同じ文字が 2 回ある: '") + cs[i] + "'");
            A.table_[c] = (short)i;
            A.chars_[i] = c;
        }
        A.size_ = (int)cs.size();
        return A;
    }
    static Alphabet ints(long long K) { return ints(0, K); }  // 整数 0〜K-1
    static Alphabet ints(long long lo, long long hi) {         // 整数 [lo, hi)
        if (hi <= lo || (unsigned long long)hi - (unsigned long long)lo > 1000000000ULL)
            str_internal::fail("Alphabet::ints", "lo < hi かつ hi - lo ≤ 10^9 で指定する（大きな値は compress_values で 0〜K-1 に）");
        Alphabet A(0);
        A.kind_ = INTS;
        A.lo_ = lo;
        A.size_ = (int)(hi - lo);
        return A;
    }

    int size() const { return size_; }
    bool is_chars() const { return kind_ == CHARS; }  // 文字の表で決めたものか（of・lower など）。false なら ints
    template <class T>
    bool contains(T x) const { return index_of(str_internal::code_of(x)) >= 0; }
    template <class T>
    int id(T x) const {  // 文字・値の番号（範囲外は止める）
        const long long v = str_internal::code_of(x);
        const int r = index_of(v);
        if (r < 0) out_of_range(v);
        return r;
    }
    long long value(int i) const {  // 番号 i の値（文字なら文字コード 0〜255）
        check_index(i);
        return kind_ == CHARS ? (long long)chars_[i] : lo_ + i;
    }
    char ch(int i) const {  // 番号 i の文字（ints で 0〜255 の外の値なら止める）
        const long long v = value(i);
        if (v < 0 || v > 255) str_internal::fail("Alphabet::ch", "番号 " + std::to_string(i) + " の値 " + std::to_string(v) + " は文字にできない");
        return (char)(unsigned char)v;
    }
    std::string chars() const {  // 文字を番号の順に並べたもの（ints なら 0〜255 の範囲のときだけ）
        std::string s;
        for (int i = 0; i < size_; ++i) s += ch(i);
        return s;
    }

   private:
    enum Kind : unsigned char { CHARS, INTS };
    explicit Alphabet(int) { table_.fill(-1); }
    static Alphabet code_range(int lo, int hi) {
        Alphabet A(0);
        A.kind_ = CHARS;
        for (int c = lo; c < hi; ++c) A.table_[c] = (short)(c - lo), A.chars_[c - lo] = (unsigned char)c;
        A.size_ = hi - lo;
        return A;
    }
    int index_of(long long v) const {
        if (kind_ == CHARS) return (unsigned long long)v < 256ULL ? table_[(std::size_t)v] : -1;
        const unsigned long long d = (unsigned long long)v - (unsigned long long)lo_;
        return d < (unsigned long long)size_ ? (int)d : -1;
    }
    void check_index(int i) const {
        if (i < 0 || i >= size_) str_internal::fail("Alphabet", "番号 " + std::to_string(i) + " は 0〜" + std::to_string(size_ - 1) + " の外");
    }
    [[noreturn]] void out_of_range(long long v) const {
        std::string what = "値 " + std::to_string(v) + " は";
        if (kind_ == CHARS && v >= 32 && v < 127) what = std::string("文字 '") + (char)v + "'（コード " + std::to_string(v) + "）は";
        str_internal::fail("Alphabet::id", what + "文字の種類に入っていない（種類の数 " + std::to_string(size_) + "）");
    }
    Kind kind_ = CHARS;
    int size_ = 0;
    long long lo_ = 0;
    std::array<short, 256> table_{};
    std::array<unsigned char, 256> chars_{};
};

template <class Seq>
std::vector<int> to_ids(const Seq& s0, const Alphabet& A = Alphabet::lower()) {  // 各文字の番号の列
    const auto& s = str_internal::seq(s0);
    std::vector<int> r(s.size());
    for (std::size_t i = 0; i < s.size(); ++i) r[i] = A.id(s[i]);
    return r;
}
inline std::string from_ids(const std::vector<int>& ids, const Alphabet& A = Alphabet::lower()) {  // 番号の列 → 文字列
    std::string s(ids.size(), '\0');
    for (std::size_t i = 0; i < ids.size(); ++i) s[i] = A.ch(ids[i]);
    return s;
}

// ---- 数の列を 0〜K-1 に圧縮（大小の順を保つ）：Alphabet::ints(K) と組み合わせる ----
struct CompressedValues {
    std::vector<int> ids;            // 各要素の番号（0〜K-1。値の小さい順）
    std::vector<long long> values;   // 番号 → 元の値（昇順。K = values.size()）
    int size() const { return (int)values.size(); }
    Alphabet alphabet() const { return Alphabet::ints(values.empty() ? 1 : (long long)values.size()); }
};
template <class Seq>
CompressedValues compress_values(const Seq& s0) {
    const auto& s = str_internal::seq(s0);
    const std::size_t n = s.size();
    std::vector<std::pair<long long, int>> vi(n);
    for (std::size_t i = 0; i < n; ++i) vi[i] = {str_internal::code_of(s[i]), (int)i};
    std::sort(vi.begin(), vi.end());
    CompressedValues c;
    c.ids.resize(n);
    for (std::size_t i = 0; i < n; ++i) {
        if (i == 0 || vi[i].first != vi[i - 1].first) c.values.push_back(vi[i].first);
        c.ids[vi[i].second] = (int)c.values.size() - 1;
    }
    return c;
}

// ---- 文字の数 ----
template <class Seq>
std::vector<int> count_chars(const Seq& s0, const Alphabet& A = Alphabet::lower()) {  // 番号ごとの個数（長さ A.size()）
    const auto& s = str_internal::seq(s0);
    std::vector<int> c(A.size(), 0);
    for (std::size_t i = 0; i < s.size(); ++i) ++c[A.id(s[i])];
    return c;
}

// ---- 区間の文字の数・次 / 前の出現を O(1) で（前処理 O(Nσ)・メモリ (N+1)σ 個の int。N = 2×10^5, σ = 26 で約 21MB） ----
//   メモリを抑えたいとき・σ が大きいときは CharPositions（O(log N)・メモリ O(N + σ)）
struct CharPrefix {
    CharPrefix() = default;
    template <class Seq>
    explicit CharPrefix(const Seq& s0, const Alphabet& A = Alphabet::lower()) : A_(A) {
        const auto& s = str_internal::seq(s0);
        n_ = (int)s.size(), sg_ = A_.size();
        if ((long long)(n_ + 1) * sg_ > 600000000LL) str_internal::fail("CharPrefix", "(N + 1)σ が 6×10^8 を超える（メモリが足りない。CharPositions を使う）");
        cnt_.assign((std::size_t)(n_ + 1) * sg_, 0);
        std::vector<int> id(n_);
        start_.assign(sg_ + 1, 0);
        for (int i = 0; i < n_; ++i) {
            const int* src = cnt_.data() + (std::size_t)i * sg_;
            int* dst = cnt_.data() + (std::size_t)(i + 1) * sg_;
            std::copy(src, src + sg_, dst);
            id[i] = A_.id(s[i]);
            ++dst[id[i]], ++start_[id[i] + 1];
        }
        for (int k = 0; k < sg_; ++k) start_[k + 1] += start_[k];
        pos_.resize(n_);
        std::vector<int> cur(start_.begin(), start_.end() - 1);
        for (int i = 0; i < n_; ++i) pos_[cur[id[i]]++] = i;
    }
    int size() const { return n_; }
    const Alphabet& alphabet() const { return A_; }
    template <class C>
    int count(int l, int r, C c) const {  // [l, r) にある c の数  O(1)
        str_internal::check_range("CharPrefix::count", l, r, n_);
        const int k = A_.id(c);
        return cnt_[(std::size_t)r * sg_ + k] - cnt_[(std::size_t)l * sg_ + k];
    }
    template <class C>
    bool has(int l, int r, C c) const { return count(l, r, c) > 0; }  // [l, r) に c があるか  O(1)
    std::vector<int> counts(int l, int r) const {  // [l, r) の番号ごとの数  O(σ)
        str_internal::check_range("CharPrefix::counts", l, r, n_);
        std::vector<int> c(sg_);
        for (int k = 0; k < sg_; ++k) c[k] = cnt_[(std::size_t)r * sg_ + k] - cnt_[(std::size_t)l * sg_ + k];
        return c;
    }
    int distinct(int l, int r) const {  // [l, r) にある文字の種類の数  O(σ)
        str_internal::check_range("CharPrefix::distinct", l, r, n_);
        int d = 0;
        for (int k = 0; k < sg_; ++k) d += cnt_[(std::size_t)r * sg_ + k] != cnt_[(std::size_t)l * sg_ + k];
        return d;
    }
    template <class C>
    int next(int i, C c) const {  // i 以上で最初の c の位置（なければ -1。0 ≤ i ≤ N）  O(1)
        check_pos("CharPrefix::next", i);
        const int k = A_.id(c), m = cnt_[(std::size_t)i * sg_ + k];
        return start_[k] + m < start_[k + 1] ? pos_[start_[k] + m] : -1;
    }
    template <class C>
    int prev(int i, C c) const {  // i より前（[0, i)）で最後の c の位置（なければ -1。0 ≤ i ≤ N）  O(1)
        check_pos("CharPrefix::prev", i);
        const int k = A_.id(c), m = cnt_[(std::size_t)i * sg_ + k];
        return m > 0 ? pos_[start_[k] + m - 1] : -1;
    }
    template <class C>
    int kth(C c, long long k) const {  // c の k 番目（0 から）の位置（なければ -1）  O(1)
        const int id = A_.id(c);
        return (k >= 0 && k < start_[id + 1] - start_[id]) ? pos_[start_[id] + k] : -1;
    }
    template <class C>
    int total(C c) const {  // 全体での c の数  O(1)
        const int k = A_.id(c);
        return start_[k + 1] - start_[k];
    }

   private:
    void check_pos(const char* func, int i) const {
        if (i < 0 || i > n_) str_internal::fail(func, "位置 " + std::to_string(i) + " は 0〜N（= " + std::to_string(n_) + "）の外");
    }
    Alphabet A_;
    int n_ = 0, sg_ = 0;
    std::vector<int> cnt_;           // cnt_[i * σ + k] = [0, i) にある番号 k の数
    std::vector<int> start_, pos_;   // 番号ごとの出現位置（pos_[start_[k] .. start_[k+1])）
};

// ---- 文字の出現位置（前処理 O(N + σ)・メモリ O(N + σ)。σ が大きくても使える） ----
struct CharPositions {
    CharPositions() = default;
    template <class Seq>
    explicit CharPositions(const Seq& s0, const Alphabet& A = Alphabet::lower()) : A_(A) {
        const auto& s = str_internal::seq(s0);
        n_ = (int)s.size();
        const int sg = A_.size();
        std::vector<int> id(n_);
        start_.assign(sg + 1, 0);
        for (int i = 0; i < n_; ++i) id[i] = A_.id(s[i]), ++start_[id[i] + 1];
        for (int k = 0; k < sg; ++k) start_[k + 1] += start_[k];
        pos_.resize(n_), occ_.resize(n_);
        std::vector<int> cur(start_.begin(), start_.end() - 1);
        for (int i = 0; i < n_; ++i) occ_[i] = cur[id[i]] - start_[id[i]], pos_[cur[id[i]]++] = i;
    }
    int size() const { return n_; }
    const Alphabet& alphabet() const { return A_; }
    template <class C>
    int next(int i, C c) const {  // i 以上で最初の c の位置（なければ -1。0 ≤ i ≤ N）  O(log N)
        check_pos("CharPositions::next", i);
        const int k = A_.id(c), j = str_internal::lower_index(pos_.data(), start_[k], start_[k + 1], i);
        return j < start_[k + 1] ? pos_[j] : -1;
    }
    template <class C>
    int prev(int i, C c) const {  // i より前（[0, i)）で最後の c の位置（なければ -1。0 ≤ i ≤ N）  O(log N)
        check_pos("CharPositions::prev", i);
        const int k = A_.id(c), j = str_internal::lower_index(pos_.data(), start_[k], start_[k + 1], i);
        return j > start_[k] ? pos_[j - 1] : -1;
    }
    template <class C>
    int kth(C c, long long k) const {  // c の k 番目（0 から）の位置（なければ -1）  O(1)
        const int id = A_.id(c);
        return (k >= 0 && k < start_[id + 1] - start_[id]) ? pos_[start_[id] + k] : -1;
    }
    template <class C>
    int count(int l, int r, C c) const {  // [l, r) にある c の数  O(log N)
        str_internal::check_range("CharPositions::count", l, r, n_);
        const int k = A_.id(c);
        return str_internal::lower_index(pos_.data(), start_[k], start_[k + 1], r) - str_internal::lower_index(pos_.data(), start_[k], start_[k + 1], l);
    }
    template <class C>
    bool has(int l, int r, C c) const { return count(l, r, c) > 0; }
    template <class C>
    int total(C c) const {  // 全体での c の数  O(1)
        const int k = A_.id(c);
        return start_[k + 1] - start_[k];
    }
    template <class C>
    std::vector<int> positions(C c) const {  // c の位置の一覧（昇順）
        const int k = A_.id(c);
        return std::vector<int>(pos_.begin() + start_[k], pos_.begin() + start_[k + 1]);
    }
    int occurrence_index(int i) const {  // s[i] が「その文字の何番目（0 から）の出現」か  O(1)
        if (i < 0 || i >= n_) str_internal::fail("CharPositions::occurrence_index", "位置 " + std::to_string(i) + " は 0〜N-1 の外");
        return occ_[i];
    }

   private:
    void check_pos(const char* func, int i) const {
        if (i < 0 || i > n_) str_internal::fail(func, "位置 " + std::to_string(i) + " は 0〜N（= " + std::to_string(n_) + "）の外");
    }
    Alphabet A_;
    int n_ = 0;
    std::vector<int> start_, pos_, occ_;
};

// ---- ランレングス圧縮：連続する同じ値を {value, len, start} にまとめる ----
template <class T>
struct RunBlock {
    T value{};
    long long len = 0;    // 続く長さ
    long long start = 0;  // 元の列での開始位置
    friend std::ostream& operator<<(std::ostream& os, const RunBlock& b) { return os << "{" << b.value << ", " << b.len << ", " << b.start << "}"; }  // デバッグ用
};
template <class Seq>
std::vector<RunBlock<str_internal::elem_t<Seq>>> run_length(const Seq& s0) {
    const auto& s = str_internal::seq(s0);
    std::vector<RunBlock<str_internal::elem_t<Seq>>> r;
    std::size_t cnt = s.empty() ? 0 : 1;
    for (std::size_t i = 1; i < s.size(); ++i) cnt += !(s[i] == s[i - 1]);
    r.reserve(cnt);
    for (std::size_t i = 0, n = s.size(); i < n;) {
        std::size_t j = i + 1;
        while (j < n && s[j] == s[i]) ++j;
        r.push_back({s[i], (long long)(j - i), (long long)i});
        i = j;
    }
    return r;
}
template <class T>
auto run_length_decode(const std::vector<RunBlock<T>>& runs) {  // 元に戻す（T = char なら string、それ以外は vector<T>）
    using R = std::conditional_t<std::is_same<T, char>::value, std::string, std::vector<T>>;
    long long total = 0;
    for (const auto& b : runs) {
        if (b.len < 0) str_internal::fail("run_length_decode", "len が負");
        total += b.len;
        if (total > 2000000000LL) str_internal::fail("run_length_decode", "元の長さが 2×10^9 を超える");
    }
    R r;
    r.reserve((std::size_t)total);
    for (const auto& b : runs) r.insert(r.end(), (std::size_t)b.len, b.value);
    return r;
}

// ---- 変換 ----
inline std::string to_lower(std::string s) {  // 英大文字を小文字に（ほかはそのまま）
    for (char& c : s)
        if ('A' <= c && c <= 'Z') c = (char)(c - 'A' + 'a');
    return s;
}
inline std::string to_upper(std::string s) {
    for (char& c : s)
        if ('a' <= c && c <= 'z') c = (char)(c - 'a' + 'A');
    return s;
}
inline std::string swap_case(std::string s) {  // 大文字と小文字を入れ替える
    for (char& c : s) {
        if ('A' <= c && c <= 'Z') c = (char)(c - 'A' + 'a');
        else if ('a' <= c && c <= 'z') c = (char)(c - 'a' + 'A');
    }
    return s;
}
template <class Seq>
str_internal::seq_t<Seq> caesar_shift(const Seq& s0, long long k, const Alphabet& A = Alphabet::lower()) {  // A の文字を番号で k ずらす（負も。A にない文字はそのまま）
    auto s = str_internal::seq_t<Seq>(str_internal::seq(s0));
    const long long sg = A.size(), sh = ((k % sg) + sg) % sg;
    for (auto& c : s)
        if (A.contains(c)) c = (std::decay_t<decltype(c)>)A.value((int)((A.id(c) + sh) % sg));
    return s;
}
// 文字の置き換えを重ねる：replace(c, d) で「今 c になっている文字を全部 d に」。apply(s) で s に当てはめる（ABC342 C）
struct CharMap {
    explicit CharMap(const Alphabet& A = Alphabet::lower()) : A_(A), f_(A.size()) {
        for (int i = 0; i < A_.size(); ++i) f_[i] = i;
    }
    template <class C, class D>
    void replace(C c, D d) {  // 今 c になっているものを全部 d に  O(σ)
        const int x = A_.id(c), y = A_.id(d);
        if (x == y) return;
        for (int& v : f_)
            if (v == x) v = y;
    }
    template <class C, class D>
    void set(C c, D d) { f_[A_.id(c)] = A_.id(d); }  // 元の文字 c の行き先を d にする（ほかは変えない）
    template <class C>
    auto operator()(C c) const {  // 元の文字 c の行き先（A にない文字はそのまま）
        return A_.contains(c) ? (C)A_.value(f_[A_.id(c)]) : c;
    }
    template <class Seq>
    str_internal::seq_t<Seq> apply(const Seq& s0) const {  // s の各文字を行き先に（A にない文字はそのまま）
        auto s = str_internal::seq_t<Seq>(str_internal::seq(s0));
        for (auto& c : s) c = (*this)(c);
        return s;
    }
    CharMap compose(const CharMap& g) const {  // 「この置き換えのあとに g」を 1 つにしたもの
        CharMap h(A_);
        for (int i = 0; i < A_.size(); ++i) h.f_[i] = g.f_[f_[i]];
        return h;
    }
    const std::vector<int>& table() const { return f_; }  // 番号 → 行き先の番号

   private:
    Alphabet A_;
    std::vector<int> f_;
};
template <class Seq>
str_internal::seq_t<Seq> reversed(const Seq& s0) {  // 反転したコピー
    const auto& s = str_internal::seq(s0);
    return str_internal::seq_t<Seq>(s.rbegin(), s.rend());
}
template <class Seq>
str_internal::seq_t<Seq> min_with_reverse(const Seq& s0) {  // s と反転のうち辞書順で小さい方（向きを区別しないときの代表。ABC310 C）
    auto s = str_internal::seq_t<Seq>(str_internal::seq(s0));
    auto t = reversed(s);
    return str_internal::lex_compare(t, s) < 0 ? t : s;
}
template <class Seq>
str_internal::seq_t<Seq> sorted_chars(const Seq& s0) {  // 並べ替えたもの（並べ替えて同じかの比較用）
    auto s = str_internal::seq_t<Seq>(str_internal::seq(s0));
    std::sort(s.begin(), s.end(), [](const auto& a, const auto& b) { return str_internal::code_of(a) < str_internal::code_of(b); });
    return s;
}
template <class Seq>
str_internal::seq_t<Seq> repeat_string(const Seq& s0, long long k) {  // s を k 回つなげる
    const auto& s = str_internal::seq(s0);
    if (k < 0) str_internal::fail("repeat_string", "k ≥ 0 で使う");
    if (!s.empty() && k > 2000000000LL / (long long)s.size()) str_internal::fail("repeat_string", "長さが 2×10^9 を超える");
    str_internal::seq_t<Seq> r;
    r.reserve((std::size_t)(k * (long long)s.size()));
    for (long long i = 0; i < k; ++i) r.insert(r.end(), s.begin(), s.end());
    return r;
}

// ---- 比較 ----
template <class A, class B>
int common_prefix_length(const A& a0, const B& b0) {  // 共通の接頭辞の長さ
    const auto& a = str_internal::seq(a0);
    const auto& b = str_internal::seq(b0);
    const std::size_t k = a.size() < b.size() ? a.size() : b.size();
    std::size_t i = 0;
    while (i < k && str_internal::code_of(a[i]) == str_internal::code_of(b[i])) ++i;
    return (int)i;
}
template <class A, class B>
int common_suffix_length(const A& a0, const B& b0) {  // 共通の接尾辞の長さ
    const auto& a = str_internal::seq(a0);
    const auto& b = str_internal::seq(b0);
    const std::size_t k = a.size() < b.size() ? a.size() : b.size();
    std::size_t i = 0;
    while (i < k && str_internal::code_of(a[a.size() - 1 - i]) == str_internal::code_of(b[b.size() - 1 - i])) ++i;
    return (int)i;
}
template <class A, class B>
int hamming_distance(const A& a0, const B& b0) {  // 同じ位置で違う文字の数（長さが同じこと）
    const auto& a = str_internal::seq(a0);
    const auto& b = str_internal::seq(b0);
    if (a.size() != b.size()) str_internal::fail("hamming_distance", "長さが違う");
    int d = 0;
    for (std::size_t i = 0; i < a.size(); ++i) d += str_internal::code_of(a[i]) != str_internal::code_of(b[i]);
    return d;
}
template <class T, class S>
bool is_subsequence(const T& t0, const S& s0) {  // t が s の部分列（飛び飛びでよい）か  O(|s|)。何度も聞くなら 06_subsequence の NextTable
    const auto& t = str_internal::seq(t0);
    const auto& s = str_internal::seq(s0);
    std::size_t j = 0;
    for (std::size_t i = 0; i < s.size() && j < t.size(); ++i)
        if (str_internal::code_of(s[i]) == str_internal::code_of(t[j])) ++j;
    return j == t.size();
}
template <class Seq>
bool is_palindrome(const Seq& s0) {  // 全体が回文か  O(N)。区間を何度も聞くなら 04_palindrome の Palindromes
    const auto& s = str_internal::seq(s0);
    for (std::size_t i = 0, n = s.size(); i < n / 2; ++i)
        if (!(s[i] == s[n - 1 - i])) return false;
    return true;
}

// ---- 区切り・つなげる ----
inline std::vector<std::string> split(const std::string& s, char sep, bool skip_empty = false) {  // sep で区切る（空の部分も残す。skip_empty で除く）
    std::vector<std::string> r;
    std::size_t b = 0;
    for (std::size_t i = 0; i <= s.size(); ++i)
        if (i == s.size() || s[i] == sep) {
            if (!(skip_empty && i == b)) r.push_back(s.substr(b, i - b));
            b = i + 1;
        }
    return r;
}
template <class V>
std::string join(const std::vector<V>& vs, const std::string& sep = " ") {  // 区切りを挟んでつなげる（数も文字列にする。答えの出力に）
    std::string r;
    for (std::size_t i = 0; i < vs.size(); ++i) {
        if (i) r += sep;
        if constexpr (std::is_same<V, std::string>::value) r += vs[i];
        else if constexpr (std::is_same<V, char>::value) r += vs[i];
        else if constexpr (std::is_integral<V>::value) r += std::to_string(vs[i]);
        else {
            std::ostringstream os;
            os << vs[i];
            r += os.str();
        }
    }
    return r;
}
#endif  // GG_STR_BASE
#ifndef GG_REVERSIBLE_DEQUE
#define GG_REVERSIBLE_DEQUE
// ---- 前後に足せて O(1) で反転できる列・回転と反転が O(1) の配列（文字列第二世代） ----
// T = char で文字列として使う（to_string・substr が string を返す）。位置は 0 始まり

// ---- ReversibleDeque<T>：前後への追加・削除と反転（環状バッファ ＋ 向きの印） ----
template <class T>
struct ReversibleDeque {
    using seq_type = std::conditional_t<std::is_same<T, char>::value, std::string, std::vector<T>>;
    ReversibleDeque() = default;
    template <class Seq, class = decltype(str_internal::seq(std::declval<const Seq&>()).size())>
    explicit ReversibleDeque(const Seq& s) { push_back(s); }

    int size() const { return n_; }
    bool empty() const { return n_ == 0; }
    void clear() { head_ = 0, n_ = 0, rev_ = false; }
    void reverse() { rev_ = !rev_; }  // 全体を反転  O(1)
    bool reversed() const { return rev_; }

    void push_back(const T& x) {  // 末尾に 1 つ  O(1)（ならし）
        T v(x);  // x がこの列の要素でも、置き場を広げる前に写しておく（d.push_back(d.front()) など）
        grow(1);
        if (!rev_) buf_[(head_ + n_) & mask_] = std::move(v);  // 置き場の後ろ
        else head_ = (head_ - 1) & mask_, buf_[head_] = std::move(v);  // 反転中は置き場の前
        ++n_;
    }
    void push_front(const T& x) {  // 先頭に 1 つ  O(1)（ならし）
        T v(x);
        grow(1);
        if (rev_) buf_[(head_ + n_) & mask_] = std::move(v);
        else head_ = (head_ - 1) & mask_, buf_[head_] = std::move(v);
        ++n_;
    }
    template <class Seq, class = decltype(str_internal::seq(std::declval<const Seq&>()).size())>
    void push_back(const Seq& s0) {  // 末尾に列をそのままの順でつなげる（自分自身も可）  O(|s|)
        if constexpr (std::is_same<Seq, ReversibleDeque>::value)
            if (&s0 == this) return push_back(to_vector());
        const auto& s = str_internal::seq(s0);
        const std::size_t m = s.size();
        grow((int)m);
        for (std::size_t i = 0; i < m; ++i) push_back((T)s[i]);
    }
    template <class Seq, class = decltype(str_internal::seq(std::declval<const Seq&>()).size())>
    void push_front(const Seq& s0) {  // 先頭に列をそのままの順でつなげる（s + 今の列。自分自身も可）  O(|s|)
        if constexpr (std::is_same<Seq, ReversibleDeque>::value)
            if (&s0 == this) return push_front(to_vector());
        const auto& s = str_internal::seq(s0);
        const std::size_t m = s.size();
        grow((int)m);
        for (std::size_t i = m; i-- > 0;) push_front((T)s[i]);
    }
    T pop_back() {  // 末尾を取り出す  O(1)
        if (n_ == 0) str_internal::fail("ReversibleDeque::pop_back", "空");
        T x = back();
        if (!rev_) --n_;
        else head_ = (head_ + 1) & mask_, --n_;
        return x;
    }
    T pop_front() {  // 先頭を取り出す  O(1)
        if (n_ == 0) str_internal::fail("ReversibleDeque::pop_front", "空");
        T x = front();
        if (rev_) --n_;
        else head_ = (head_ + 1) & mask_, --n_;
        return x;
    }
    const T& front() const { return at(0, "ReversibleDeque::front"); }
    const T& back() const { return at(n_ - 1, "ReversibleDeque::back"); }
    const T& operator[](int i) const { return at(i, "ReversibleDeque::operator[]"); }  // 先頭から i 番目  O(1)
    T& operator[](int i) { return const_cast<T&>(at(i, "ReversibleDeque::operator[]")); }
    std::vector<T> to_vector() const {  // 先頭から順に  O(N)
        std::vector<T> r(n_);
        for (int i = 0; i < n_; ++i) r[i] = buf_[phys(i)];
        return r;
    }
    std::string to_string() const {  // 文字列で（T = char）  O(N)
        static_assert(std::is_same<T, char>::value, "to_string は T = char のときだけ（ほかは to_vector）");
        std::string r(n_, '\0');
        for (int i = 0; i < n_; ++i) r[i] = buf_[phys(i)];
        return r;
    }
    seq_type substr(int pos, int len) const {  // [pos, pos + len) を取り出す（T = char なら string）  O(len)
        str_internal::check_range("ReversibleDeque::substr", pos, (long long)pos + len, n_);
        seq_type r(len, T());
        for (int i = 0; i < len; ++i) r[i] = buf_[phys(pos + i)];
        return r;
    }

   private:
    int phys(int i) const { return (head_ + (rev_ ? n_ - 1 - i : i)) & mask_; }  // 先頭から i 番目の場所
    const T& at(int i, const char* func) const {
        if (i < 0 || i >= n_) str_internal::fail(func, "位置 " + std::to_string(i) + " は 0〜" + std::to_string(n_ - 1) + " の外");
        return buf_[phys(i)];
    }
    void grow(int add) {
        if ((long long)n_ + add <= (long long)buf_.size()) return;
        std::size_t cap = buf_.empty() ? 16 : buf_.size();
        while (cap < (std::size_t)n_ + (std::size_t)add) cap <<= 1;
        std::vector<T> nb(cap);
        for (int i = 0; i < n_; ++i) nb[i] = buf_[(head_ + i) & mask_];  // 置き場の順のまま詰める（向きの印はそのまま）
        buf_.swap(nb);
        head_ = 0, mask_ = (int)cap - 1;
    }
    std::vector<T> buf_;
    int head_ = 0, n_ = 0, mask_ = -1;
    bool rev_ = false;
};

// ---- RotatableArray<T>：大きさが決まった配列の回転と反転を O(1)（ずれと向きを持つ） ----
template <class T>
struct RotatableArray {
    using seq_type = std::conditional_t<std::is_same<T, char>::value, std::string, std::vector<T>>;
    RotatableArray() = default;
    template <class Seq>
    explicit RotatableArray(const Seq& s0) {
        const auto& s = str_internal::seq(s0);
        a_.assign(s.begin(), s.end());
        n_ = (long long)a_.size();
    }
    int size() const { return (int)n_; }
    void rotate_left(long long k) {  // 左に k 回す（新しい [i] = 元の [i + k]。負なら右）  O(1)
        if (n_ == 0) return;
        k %= n_;
        off_ = ((off_ + dir_ * k) % n_ + n_) % n_;
    }
    void rotate_right(long long k) { rotate_left(-(k % (n_ == 0 ? 1 : n_))); }  // 右に k 回す（新しい [i] = 元の [i - k]）  O(1)
    void reverse() {  // 全体を反転  O(1)
        if (n_ == 0) return;
        off_ = ((off_ + dir_ * (n_ - 1)) % n_ + n_) % n_;
        dir_ = -dir_;
    }
    const T& operator[](int i) const { return a_[phys(i, "RotatableArray::operator[]")]; }  // 今の i 番目  O(1)
    T& operator[](int i) { return a_[phys(i, "RotatableArray::operator[]")]; }
    const T& get(int i) const { return a_[phys(i, "RotatableArray::get")]; }
    void set(int i, const T& v) { a_[phys(i, "RotatableArray::set")] = v; }
    std::vector<T> to_vector() const {  // 今の順に  O(N)
        std::vector<T> r(n_);
        for (int i = 0; i < (int)n_; ++i) r[i] = (*this)[i];
        return r;
    }
    std::string to_string() const {  // 文字列で（T = char）  O(N)
        static_assert(std::is_same<T, char>::value, "to_string は T = char のときだけ（ほかは to_vector）");
        std::string r((std::size_t)n_, '\0');
        for (int i = 0; i < (int)n_; ++i) r[i] = (*this)[i];
        return r;
    }

   private:
    std::size_t phys(int i, const char* func) const {
        if (i < 0 || i >= n_) str_internal::fail(func, "位置 " + std::to_string(i) + " は 0〜" + std::to_string(n_ - 1) + " の外");
        long long p = off_ + dir_ * (long long)i;
        if (p >= n_) p -= n_;
        else if (p < 0) p += n_;
        return (std::size_t)p;
    }
    std::vector<T> a_;
    long long n_ = 0, off_ = 0, dir_ = 1;
};
#endif  // GG_REVERSIBLE_DEQUE

namespace old_rh {
struct RollingHash {
    using ull = unsigned long long;
    static constexpr ull mod = (1ULL << 61) - 1;

    /**
     * @brief ハッシュ情報
     *
     * value : ハッシュ値
     * left  : 元文字列での開始位置（0-index）
     * right : 元文字列での終了位置（0-index, 半開区間）
     * len   : 部分文字列長
     *
     * ※ concat() で生成した Hash は len のみ有効で、
     *    left==-1 の場合は元文字列に対応しない Hash を表す。
     */
    struct Hash {
        ull value;
        int left;
        int right;
        int len;

        operator ull() const { return value; }

        bool operator==(const Hash& rhs) const {
            return value == rhs.value;
        }

        bool operator!=(const Hash& rhs) const {
            return value != rhs.value;
        }

        bool operator<(const Hash& rhs) const {
            return value < rhs.value;
        }
    };

    static ull base;
    static bool initialized;
    static vector<ull> power;

    vector<ull> h;

    /**
     * @brief base をランダムに初期化
     *
     * 初回のみ呼ばれ、base と power を初期化する。
     * コンストラクタで呼ばれるため個別での実行不要
     *
     * Time Complexity:
     *   O(1)
     */
    static void init() {
        if (initialized) return;
        initialized = true;

        mt19937_64 rng(
            chrono::steady_clock::now().time_since_epoch().count()
        );
        base = uniform_int_distribution<ull>(
            100000000ULL,
            1000000000ULL
        )(rng);

        power = {1};
    }

    /**
     * @brief 文字列から Rolling Hash を構築
     *
     * @param s 元文字列
     *
     * Time Complexity:
     *   O(N)
     */
    RollingHash(const string& s) {
        init();

        int n = (int)s.size();
        ensure_power(n);

        h.assign(n + 1, 0);

        for (int i = 0; i < n; i++) {
            h[i + 1] = add(
                mul(h[i], base),
                (ull)(unsigned char)s[i] + 1
            );
        }
    }

    /**
     * @brief 文字列長を取得
     *
     * @return 文字列長
     *
     * Time Complexity:
     *   O(1)
     */
    int size() const {
        return (int)h.size() - 1;
    }

    /**
     * @brief 部分文字列 [l,r) のハッシュを取得
     *
     * @param l 開始位置（0-index）
     * @param r 終了位置（0-index, 半開区間）
     * @return Hash
     *
     * Time Complexity:
     *   O(1)
     */
    Hash get(int l, int r) const {
        return {
            add(h[r], mod - mul(h[l], power[r - l])),
            l,
            r,
            r - l
        };
    }

    /**
     * @brief 長さ指定で部分文字列のハッシュを取得
     *
     * s[l,l+len) のハッシュを返す。
     *
     * @param l 開始位置
     * @param len 長さ
     * @return Hash
     *
     * Time Complexity:
     *   O(1)
     */
    Hash get_len(int l, int len) const {
        return get(l, l + len);
    }

    /**
     * @brief 接尾辞のハッシュを取得
     *
     * s[l,n) のハッシュを返す。
     *
     * @param l 開始位置
     * @return Hash
     *
     * Time Complexity:
     *   O(1)
     */
    Hash get_from(int l = 0) const {
        return get(l, size());
    }

    /**
     * @brief 空文字列のハッシュを取得
     *
     * concat() の単位元として利用できる。
     *
     * @return Hash
     *
     * Time Complexity:
     *   O(1)
     */
    static Hash empty_hash() {
        return {
            0,
            -1,
            -1,
            0
        };
    }

    /**
     * @brief 2つのハッシュを連結
     *
     * a + b のハッシュを返す。
     * 異なる RollingHash 同士でも利用可能。
     *
     * ※返却される Hash は len のみ有効で、
     *   left==-1 のため元文字列には対応しない。
     *
     * @param a 前半
     * @param b 後半
     * @return Hash
     *
     * Time Complexity:
     *   O(1)
     */
    static Hash concat(const Hash& a, const Hash& b) {
        ensure_power(b.len);

        return {
            add(
                mul(a.value, power[b.len]),
                b.value
            ),
            -1,
            -1,
            a.len + b.len
        };
    }

    /**
     * @brief 最長共通接頭辞長（LCP）を取得
     *
     * a[i,n) と b[j,n) の LCP を返す。
     *
     * @param a 文字列A
     * @param i 開始位置
     * @param b 文字列B
     * @param j 開始位置
     * @return LCP長
     *
     * Time Complexity:
     *   O(log N)
     */
    static int lcp(
        const RollingHash& a, int i,
        const RollingHash& b, int j
    ) {
        int low = 0;
        int high = min(a.size() - i, b.size() - j) + 1;

        while (high - low > 1) {
            int mid = (low + high) / 2;

            if (a.get(i, i + mid) == b.get(j, j + mid))
                low = mid;
            else
                high = mid;
        }

        return low;
    }

    /**
     * @brief 接尾辞同士を辞書順比較
     *
     * a[i,n) < b[j,n) なら true を返す。
     *
     * @param a 文字列A
     * @param i 開始位置
     * @param b 文字列B
     * @param j 開始位置
     * @return 比較結果
     *
     * Time Complexity:
     *   O(log N)
     */
    static bool less_than(
        const RollingHash& a, int i,
        const RollingHash& b, int j
    ) {
        int l = lcp(a, i, b, j);

        if (i + l == a.size() || j + l == b.size()) {
            return (a.size() - i) < (b.size() - j);
        }

        return a.get(i + l, i + l + 1).value
             < b.get(j + l, j + l + 1).value;
    }

    #ifdef LOCAL
    /**
     * @brief Hash の情報を表示（デバッグ用）
     *
     * ハッシュ値・区間・長さ・対応する部分文字列を表示する。
     * 提出時（-DLOCAL なし）は何もしない。
     *
     * @param s 元文字列
     * @param hash 表示する Hash
     *
     * Time Complexity:
     *   O(len)
     */
    void print_hash(const string& s, const Hash& hash) const {
        cout << DBG_COLOR;
        cout << "hash  = " << hash.value << '\n';
        cout << "left  = " << hash.left << '\n';
        cout << "right = " << hash.right << '\n';
        cout << "len   = " << hash.len << '\n';

        if (hash.left == -1) {
            cout << "str   = (concat 等により元文字列に対応しません)\n";
        } else {
            cout << "str   = "
                 << s.substr(hash.left, hash.len) << '\n';
        }

        cout << DBG_RESET;
    }
    #else
    /**
     * @brief Hash の情報を表示（デバッグ用）
     *
     * ハッシュ値・区間・長さ・対応する部分文字列を表示する。
     * 提出時（-DLOCAL なし）は何もしない。
     *
     * @param s 元文字列
     * @param hash 表示する Hash
     *
     * Time Complexity:
     *   O(len)
     */
    void print_hash(const string&, const Hash&) const {}
    #endif

    static ull mul(ull a, ull b) {
        __uint128_t t = (__uint128_t)a * b;
        t = (t >> 61) + (t & mod);
        if (t >= mod) t -= mod;
        return (ull)t;
    }

    static ull add(ull a, ull b) {
        a += b;
        if (a >= mod) a -= mod;
        return a;
    }

    static void ensure_power(int n) {
        while ((int)power.size() <= n) {
            power.push_back(
                mul(power.back(), base)
            );
        }
    }
};

RollingHash::ull RollingHash::base;
bool RollingHash::initialized = false;
vector<RollingHash::ull> RollingHash::power;
}
namespace old_man {
pair<vector<int>,vector<int>> manacher_raw(const string& S){
    int N=sizei(S);

    vector<int> odd(N),even(N);

    // odd
    for(int l=0,r=-1,i=0;i<N;i++){
        int k=(i>r?1:min(odd[l+r-i],r-i+1));
        while(0<=i-k&&i+k<N&&S[i-k]==S[i+k])k++;
        odd[i]=k--;
        if(i+k>r){
            l=i-k;
            r=i+k;
        }
    }

    // even
    for(int l=0,r=-1,i=0;i<N;i++){
        int k=(i>r?0:min(even[l+r-i+1],r-i+1));
        while(0<=i-k-1&&i+k<N&&S[i-k-1]==S[i+k])k++;
        even[i]=k--;
        if(i+k>r){
            l=i-k-1;
            r=i+k;
        }
    }

    return {odd,even};
}

// Manacher Algorithm
//
// 構築 O(N)
//
// odd_radius[i]  : iを中心とする奇数長回文の半径
// even_radius[i] : i-1,i間を中心とする偶数長回文の半径
// odd_length[i]  : 奇数長回文の最大長
// even_length[i] : 偶数長回文の最大長
//
// is_palindrome(l,r) : S[l,r]が回文か判定 O(1)
struct Manacher{
    vector<int> odd_radius;
    vector<int> even_radius;

    vector<int> odd_length;
    vector<int> even_length;

    Manacher(const string& S){
        tie(odd_radius,even_radius)=manacher_raw(S);

        int N=sizei(S);
        odd_length.resize(N);
        even_length.resize(N);

        rep(i,N){
            odd_length[i]=odd_radius[i]*2-1;
            even_length[i]=even_radius[i]*2;
        }
    }

    // S[l,r]が回文か
    bool is_palindrome(int l,int r)const{
        int len=r-l+1;
        if(len&1){
            int c=(l+r)/2;
            return odd_radius[c]>=len/2+1;
        }else{
            int c=(l+r+1)/2;
            return even_radius[c]>=len/2;
        }
    }
};
}
namespace old_trie {
//stringの追加削除、完全一致検索、prefix検索、辞書順k番目取得、辞書順取得を高速に
// ★ 文字種を変える
const int ALPHA = 26; // ★ 'a'～'z'

struct Trie {
    struct Node {
        int nxt[ALPHA];
        long long cnt;   // このノードを通る文字列数
        long long end;   // このノードで終わる文字列数
        Node(){
            fill(nxt, nxt+ALPHA, -1);
            cnt = end = 0;
        }
    };

    vector<Node> t;

    Trie(){
        t.emplace_back(); // root
    }

    // ===== 挿入・削除 =====
    // delta = +1 : 追加,  -1 : 削除
    void insert(const string& s, int delta){
        int v = 0;
        for(char c : s){
            int id = c - 'A'; // ★
            if(t[v].nxt[id] == -1){
                if(delta == -1) return; // 無いものは消さない
                t[v].nxt[id] = t.size();
                t.emplace_back();
            }
            v = t[v].nxt[id];
            t[v].cnt += delta;
        }
        t[v].end += delta;
    }

    // ===== 完全一致 =====
    bool find(const string& s){
        int v = 0;
        for(char c : s){
            int id = c - 'A'; // ★
            if(t[v].nxt[id] == -1) return false;
            v = t[v].nxt[id];
        }
        return t[v].end > 0;
    }

    // ===== prefix存在 =====
    bool starts_with(const string& s){
        int v = 0;
        for(char c : s){
            int id = c - 'A'; // ★
            if(t[v].nxt[id] == -1) return false;
            v = t[v].nxt[id];
        }
        return true;
    }

    // ===== k番目に小さい文字列（辞書順）=====
    // ★ 1-indexed
    string kth(long long k){
        string res;
        int v = 0;
        while(true){
            // このノードで終わる文字列
            if(t[v].end){
                if(k <= t[v].end) return res;
                k -= t[v].end;
            }
            // 子ノードを辞書順に
            for(int c=0;c<ALPHA;c++){
                int to = t[v].nxt[c];
                if(to == -1) continue;
                if(t[to].cnt >= k){
                    res.push_back('A' + c); // ★
                    v = to;
                    break;
                }else{
                    k -= t[to].cnt;
                }
            }
        }
    }

    void dfs(int v, string& cur, vector<string>& res){
        // このノードで終わる文字列
        for(int i=0;i<t[v].end;i++){
            res.push_back(cur);
        }
        // 子ノード
        for(int c=0;c<ALPHA;c++){
            int to = t[v].nxt[c];
            if(to == -1) continue;
            cur.push_back('A' + c); // ★
            dfs(to, cur, res);
            cur.pop_back();
        }
    }

    // ===== 辞書順で全て取得 =====
    vector<string> get_all(){
        vector<string> res;
        string cur;
        dfs(0, cur, res);
        return res;
    }

};
}
namespace old_lis {
//最長増加部分列(LIS)を返却
vector<int> longest_increasing_subsequence(const vector<int>& a){
    int n = a.size();
    if(n==0) return {};

    vector<int> dp;                 // 値を保持する最小末尾列
    vector<int> idx;                // dp の各要素が指す a のインデックス
    vector<int> parent(n, -1);      // 復元用の親ポインタ

    for(int i=0;i<n;i++){
    int left = -1,right = (int)dp.size(); // right は「候補の最小インデックス」
    while(right - left > 1){
        int mid = (left + right) / 2;
        if(dp[mid] < a[i]) left = mid; // dp[mid] < target のとき左を動かす
        else right = mid;              // dp[mid] >= target のとき右を動かす
    }
    int pos = right;

    if(pos == (int)dp.size()){
        dp.push_back(a[i]);
        idx.push_back(i);
    } else {
        dp[pos] = a[i];
        idx[pos] = i;
    }
        if(pos > 0) parent[i] = idx[pos-1];
    }

    // 最後の要素のインデックスは idx.back()
    int k = dp.size();
    vector<int> lis(k);
    int cur = idx[k-1];
    for(int i=k-1;i>=0;i--){
        lis[i] = a[cur];
        cur = parent[cur];
    }
    return lis;
}
}
#undef rep
static double now_ms() { return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now().time_since_epoch()).count(); }
template <class F> static double best_of(int k, F f) { double b = 1e18; for (int i = 0; i < k; ++i) { double t = now_ms(); f(); b = std::min(b, now_ms() - t); } return b; }
int main() {
    std::mt19937_64 rng(7);
    long long chk = 0;
    auto mk = [&](int n, int sg) { std::string s(n, 'a'); for (auto& c : s) c = (char)('a' + rng() % sg); return s; };
    // ---- 接尾辞配列・LCP・Z（ACL と比べる）----
    for (int sg : {26, 2, 1}) {
        const std::string s = mk(1000000, sg);
        std::vector<int> sa1, sa2;
        double a = best_of(3, [&] { sa1 = atcoder::suffix_array(s); });
        double b = best_of(3, [&] { sa2 = make_suffix_array(s); });
        if (sa1 != sa2) { puts("SA MISMATCH"); return 1; }
        std::vector<int> l1, l2;
        double c = best_of(3, [&] { l1 = atcoder::lcp_array(s, sa1); });
        double d = best_of(3, [&] { l2 = make_lcp_array(s, sa2); });
        if (l1 != l2) { puts("LCP MISMATCH"); return 1; }
        printf("SA  N=1e6 sigma=%2d : ACL suffix_array %6.1f ms / make_suffix_array %6.1f ms | ACL lcp_array %6.1f ms / make_lcp_array %6.1f ms\n", sg, a, b, c, d);
    }
    {
        std::vector<long long> v(1000000);
        for (auto& x : v) x = (long long)(rng() % 2000000000) - 1000000000;
        std::vector<int> sa1, sa2;
        double a = best_of(3, [&] { sa1 = atcoder::suffix_array(v); });
        double b = best_of(3, [&] { sa2 = make_suffix_array(v); });
        if (sa1 != sa2) { puts("SA(ll) MISMATCH"); return 1; }
        printf("SA  N=1e6 long long  : ACL suffix_array %6.1f ms / make_suffix_array %6.1f ms\n", a, b);
    }
    for (int sg : {26, 2, 1}) {
        const std::string s = mk(10000000, sg);
        std::vector<int> z1, z2;
        double a = best_of(3, [&] { z1 = atcoder::z_algorithm(s); });
        double b = best_of(3, [&] { z2 = z_array(s); });
        if (z1 != z2) { puts("Z MISMATCH"); return 1; }
        printf("Z   N=1e7 sigma=%2d : ACL z_algorithm %6.1f ms / z_array %6.1f ms\n", sg, a, b);
    }
    // ---- ローリングハッシュ（旧 RollingHashTest2 と比べる）----
    {
        const std::string s = mk(10000000, 26);
        std::vector<std::pair<int, int>> q(10000000);
        for (auto& [l, r] : q) { l = (int)(rng() % 10000001); r = (int)(rng() % 10000001); if (l > r) std::swap(l, r); }
        double a = best_of(3, [&] { old_rh::RollingHash H(s); chk += H.size(); });
        double b = best_of(3, [&] { RollingHash H(s); chk += H.size(); });
        old_rh::RollingHash OH(s);
        RollingHash NH(s);
        double c = best_of(3, [&] { unsigned long long x = 0; for (auto [l, r] : q) x += OH.get(l, r).value; chk += (long long)(x & 1); });
        double d = best_of(3, [&] { unsigned long long x = 0; for (auto [l, r] : q) x += NH.get(l, r).value; chk += (long long)(x & 1); });
        const std::string t = mk(1000000, 2);
        old_rh::RollingHash OT(t);
        RollingHash NT(t);
        std::vector<std::pair<int, int>> p(1000000);
        for (auto& [i, j] : p) i = (int)(rng() % 1000000), j = (int)(rng() % 1000000);
        double e = best_of(3, [&] { long long x = 0; for (auto [i, j] : p) x += old_rh::RollingHash::lcp(OT, i, OT, j); chk += x; });
        double f = best_of(3, [&] { long long x = 0; for (auto [i, j] : p) x += NT.lcp(i, j); chk += x; });
        printf("RH  構築 N=1e7     : 旧 %6.1f ms / 新 %6.1f ms | get x1e7 : 旧 %6.1f ms / 新 %6.1f ms | lcp x1e6 (N=1e6, sigma=2) : 旧 %6.1f ms / 新 %6.1f ms\n", a, b, c, d, e, f);
    }
    // ---- 回文（旧 ManacherTemplate と比べる）----
    {
        const std::string s = mk(10000000, 2);
        double a = best_of(3, [&] { old_man::Manacher M(s); chk += M.odd_radius[5]; });
        double b = best_of(3, [&] { Palindromes P(s); chk += P.size(); });
        old_man::Manacher OM(s);
        Palindromes NP(s);
        std::vector<std::pair<int, int>> q(10000000);
        for (auto& [l, r] : q) { l = (int)(rng() % 10000000); r = l + (int)(rng() % 20); if (r > 10000000) r = 10000000; }
        double c = best_of(3, [&] { long long x = 0; for (auto [l, r] : q) x += l < r ? OM.is_palindrome(l, r - 1) : 1; chk += x; });
        double d = best_of(3, [&] { long long x = 0; for (auto [l, r] : q) x += NP.is_palindrome(l, r); chk += x; });
        long long xo = 0, xn = 0;
        for (auto [l, r] : q) xo += l < r ? OM.is_palindrome(l, r - 1) : 1, xn += NP.is_palindrome(l, r);
        if (xo != xn) { puts("PAL MISMATCH"); return 1; }
        printf("PAL 構築 N=1e7     : 旧 Manacher %6.1f ms / Palindromes %6.1f ms | is_palindrome x1e7 : 旧 %6.1f ms / 新 %6.1f ms\n", a, b, c, d);
    }
    // ---- Trie（旧 ALDS14C と比べる。英大文字）----
    {
        std::vector<std::string> vs;
        long long tot = 0;
        while (tot < 1000000) { std::string w = mk(1 + (int)(rng() % 20), 26); for (auto& ch : w) ch = (char)(ch - 'a' + 'A'); tot += (long long)w.size(); vs.push_back(w); }
        double a = best_of(3, [&] { old_trie::Trie T; for (auto& w : vs) T.insert(w, 1); long long x = 0; for (auto& w : vs) x += T.find(w); chk += x; });
        double b = best_of(3, [&] { Trie T(Alphabet::upper()); for (auto& w : vs) T.insert(w); long long x = 0; for (auto& w : vs) x += T.contains(w); chk += x; });
        old_trie::Trie OT; Trie NT(Alphabet::upper());
        for (auto& w : vs) OT.insert(w, 1), NT.insert(w);
        double c = best_of(3, [&] { long long x = 0; for (int k = 0; k < 100000; ++k) x += (long long)OT.kth(1 + (long long)(k * 7919LL % (long long)vs.size())).size(); chk += x; });
        double d = best_of(3, [&] { long long x = 0; for (int k = 0; k < 100000; ++k) x += (long long)NT.kth((long long)(k * 7919LL % (long long)vs.size())).size(); chk += x; });
        if (OT.kth(1000) != NT.kth(999)) { puts("TRIE MISMATCH"); return 1; }
        printf("TRIE 文字 1e6 insert + 検索 : 旧 %6.1f ms / 新 %6.1f ms | kth x1e5 : 旧 %6.1f ms / 新 %6.1f ms\n", a, b, c, d);
    }
    // ---- LIS（旧 A24 と比べる）----
    {
        std::vector<int> v(1000000);
        for (auto& x : v) x = (int)(rng() % 1000000000);
        std::vector<int> r1, r2;
        double a = best_of(3, [&] { r1 = old_lis::longest_increasing_subsequence(v); });
        double b = best_of(3, [&] { r2 = lis_restore(v); });
        if (r1.size() != r2.size()) { puts("LIS MISMATCH"); return 1; }
        double c = best_of(3, [&] { chk += lis_length(v); });
        printf("LIS N=1e6 復元     : 旧 A24 %6.1f ms / lis_restore %6.1f ms（lis_length %6.1f ms）\n", a, b, c);
    }
    printf("check %lld\n", chk % 7);
}
