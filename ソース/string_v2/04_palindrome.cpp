// =============================================================================
//  回文（Manacher・区間の回文判定 O(1)・最長の回文・回文の個数・回文にする・最長回文部分列・回文分割・回文部分列の数）   --- 文字列第二世代
// -----------------------------------------------------------------------------
//  分類       : 04 回文
//  依存       : 00_str_base（同梱済み）
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）
//  文字コード : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  列は string でも vector<int>・vector<long long> でもよい。範囲は [l, r)。全体が回文かだけなら 00_str_base の is_palindrome(s)
//  中心の番号 c（0 ≤ c < 2N - 1）：c = 2i は s[i] が中心（奇数長）、c = 2i - 1 は s[i-1] と s[i] の間（偶数長）
//  「部分文字列」は連続する部分、「部分列」は飛び飛びに選んだもの
//
//  機能一覧                                                                       計算量
//    center_palindrome_lengths(s)          中心 c ごとの最長の回文の長さ（Manacher の生の値）         O(N)
//    Palindromes P(s)                      前処理 O(N)
//      P.is_palindrome(l, r)               s[l, r) が回文か                                          O(1)
//      P.radius_odd(i) / radius_even(i)    s[i] が中心（長さ 2r-1）/ s[i-1] と s[i] の間（長さ 2r）の半径   O(1)
//      P.length_at_center(c) / span_at_center(c)   中心 c の最長の回文の長さ / 位置
//      P.longest()                         最長の回文の部分文字列（StrSpan。同じ長さなら左）           O(N)
//      P.count_all()                       回文である部分文字列の数（位置が違えば別）                   O(N)
//      P.longest_prefix_palindrome() / longest_suffix_palindrome()   最長の回文の接頭辞 / 接尾辞の長さ  O(N)
//      P.longest_starting_at() / longest_ending_at()   各位置から始まる / で終わる最長の回文の長さ      O(N)
//    min_palindrome_with_prefix(s) / min_palindrome_with_suffix(s)   s で始まる / 終わる最短の回文     O(N)
//    min_changes_to_palindrome(s)          書き換えて回文にする最小の回数                               O(N)
//    can_rearrange_to_palindrome(s)        並べ替えて回文にできるか                                     O(N log N)
//    longest_palindromic_subsequence(s)    最長回文部分列を 1 つ（N ≤ 10^4）                            O(N^2)
//    longest_palindromic_subsequence_length(s)   その長さ（メモリ O(N)）                                O(N^2)
//    palindrome_partition_min(s)           回文に分ける最小の個数の分け方（vector<StrSpan>）              O(N^2)
//    count_palindromic_subsequences<T>(s)  回文である部分列の数（位置が違えば別。T は mint・long long）   O(N^2)
//    （異なる回文の種類・回文分割を N log N で：04_palindromic_tree（☆）。書き換えがあるとき：09_mutable_string（☆））
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    pal_test::random_test()  … 全部の区間・全部の部分列を調べる愚直解と全機能を突き合わせるランダムテスト
//    pal_test::benchmark()    … 速度計測
//    pal_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
// =============================================================================
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

// ======================= ここからテンプレート本体（コピー範囲） =======================
#ifndef GG_STR_BASE
#define GG_STR_BASE
// ---- 文字列基盤（文字列第二世代）  文字の種類 Alphabet・文字の数と位置・ランレングス・変換・比較・区切り・隣が同じにならない並べ替え ----
// 列（string / vector<整数>）を受け取る関数はどちらでもよい。文字列リテラル（"abc"）も渡せる
// 範囲は半開区間 [l, r)（0 始まり）。位置が「ない」ときは -1。使える範囲の外の引数はメッセージを出して止める
// 数（count_chars・CharPrefix::count など）は long long で返す（受け取る変数も ll に。chmax(int&, ll) はコンパイルエラーになる）
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
inline std::string seq(char c) { return std::string(1, c); }  // 1 文字も長さ 1 の列として（find_all(s, 'a') など）
template <class T>
using seq_t = std::decay_t<decltype(seq(std::declval<const T&>()))>;
template <class T>
using elem_t = std::decay_t<decltype(std::declval<const seq_t<T>&>()[0])>;
// 要素を「値」に（char は unsigned char として 0〜255。負の添字にしない）。整数の型だけ（unsigned long long は 2^63 未満の値で使う）
template <class T>
long long code_of(T x) {
    static_assert(std::is_integral<T>::value, "文字・整数だけを使える（文字列リテラル \"c\" ではなく 'c'。ポインタ・小数は不可）");
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
// 数え上げの掛け算：T が整数型なら溢れたら止める
template <class T>
T checked_mul(const T& a, const T& b, const char* func) {
    if constexpr (std::is_integral<T>::value) {
        T r;
        if (__builtin_mul_overflow(a, b, &r)) fail(func, "値が型に収まらない（答えを mod で求めるなら math_v2 の ModInt998 などを使う）");
        return r;
    } else {
        return a * b;
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
    Alphabet with(const std::string& extra) const {  // 後ろに文字を足したもの（Alphabet::lower().with("@#") など。同じ文字があると止まる）
        if (kind_ != CHARS) str_internal::fail("Alphabet::with", "ints には足せない");
        return of(chars() + extra);
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
        str_internal::fail("Alphabet::id", what + "文字の種類に入っていない（種類の数 " + std::to_string(size_) +
                                               "。既定は英小文字。Alphabet::upper()・digits()・letters()・alnum()・bytes()・of(\"…\")・lower().with(\"…\") などを渡す）");
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
template <class I>
std::string from_ids(const std::vector<I>& ids, const Alphabet& A = Alphabet::lower()) {  // 番号の列（int・long long）→ 文字列
    std::string s(ids.size(), '\0');
    for (std::size_t i = 0; i < ids.size(); ++i) {
        const long long v = str_internal::code_of(ids[i]);
        if (v < 0 || v >= A.size()) str_internal::fail("from_ids", "番号 " + std::to_string(v) + " は 0〜σ-1（σ = " + std::to_string(A.size()) + "）の外");
        s[i] = A.ch((int)v);
    }
    return s;
}
inline std::string from_ids(const std::vector<int>& ids, const Alphabet& A = Alphabet::lower()) { return from_ids<int>(ids, A); }  // from_ids({7, 8}) も

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
std::vector<long long> count_chars(const Seq& s0, const Alphabet& A = Alphabet::lower()) {  // 番号ごとの個数（長さ A.size()。long long なので積も溢れにくい）
    const auto& s = str_internal::seq(s0);
    std::vector<long long> c(A.size(), 0);
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
        if ((long long)(n_ + 1) * sg_ > 200000000LL) str_internal::fail("CharPrefix", "(N + 1)σ が 2×10^8 を超える（約 800MB。メモリが足りない。CharPositions を使う）");
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
    long long count(int l, int r, C c) const { return count_id(l, r, A_.id(c)); }  // [l, r) にある c の数  O(1)
    long long count_id(int l, int r, int k) const {                                 // 番号 k（0〜σ-1）で  O(1)
        str_internal::check_range("CharPrefix::count", l, r, n_), check_id("CharPrefix::count_id", k);
        return cnt_[(std::size_t)r * sg_ + k] - cnt_[(std::size_t)l * sg_ + k];
    }
    template <class C>
    bool has(int l, int r, C c) const { return count(l, r, c) > 0; }  // [l, r) に c があるか  O(1)
    std::vector<long long> counts(int l, int r) const {  // [l, r) の番号ごとの数  O(σ)
        str_internal::check_range("CharPrefix::counts", l, r, n_);
        std::vector<long long> c(sg_);
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
    int next(int i, C c) const { return next_id(i, A_.id(c)); }  // i 以上で最初の c の位置（なければ -1。0 ≤ i ≤ N）  O(1)
    int next_id(int i, int k) const {
        check_pos("CharPrefix::next", i), check_id("CharPrefix::next_id", k);
        const int m = cnt_[(std::size_t)i * sg_ + k];
        return start_[k] + m < start_[k + 1] ? pos_[start_[k] + m] : -1;
    }
    template <class C>
    int prev(int i, C c) const { return prev_id(i, A_.id(c)); }  // i より前（[0, i)）で最後の c の位置（なければ -1。0 ≤ i ≤ N）  O(1)
    int prev_id(int i, int k) const {
        check_pos("CharPrefix::prev", i), check_id("CharPrefix::prev_id", k);
        const int m = cnt_[(std::size_t)i * sg_ + k];
        return m > 0 ? pos_[start_[k] + m - 1] : -1;
    }
    template <class C>
    int kth(C c, long long k) const { return kth_id(A_.id(c), k); }  // c の k 番目（0 から）の位置（なければ -1）  O(1)
    int kth_id(int id, long long k) const {
        check_id("CharPrefix::kth_id", id);
        return (k >= 0 && k < start_[id + 1] - start_[id]) ? pos_[start_[id] + k] : -1;
    }
    template <class C>
    long long total(C c) const { return total_id(A_.id(c)); }  // 全体での c の数  O(1)
    long long total_id(int k) const {                          // 番号 k で（rep(k, σ) で回すときは id 版を使う）
        check_id("CharPrefix::total_id", k);
        return start_[k + 1] - start_[k];
    }

   private:
    void check_pos(const char* func, int i) const {
        if (i < 0 || i > n_) str_internal::fail(func, "位置 " + std::to_string(i) + " は 0〜N（= " + std::to_string(n_) + "）の外");
    }
    void check_id(const char* func, int k) const {
        if (k < 0 || k >= sg_) str_internal::fail(func, "番号 " + std::to_string(k) + " は 0〜σ-1（σ = " + std::to_string(sg_) + "）の外");
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
    long long count(int l, int r, C c) const {  // [l, r) にある c の数  O(log N)
        str_internal::check_range("CharPositions::count", l, r, n_);
        const int k = A_.id(c);
        return str_internal::lower_index(pos_.data(), start_[k], start_[k + 1], r) - str_internal::lower_index(pos_.data(), start_[k], start_[k + 1], l);
    }
    template <class C>
    bool has(int l, int r, C c) const { return count(l, r, c) > 0; }
    template <class C>
    long long total(C c) const {  // 全体での c の数  O(1)
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
inline char to_lower(char c) { return ('A' <= c && c <= 'Z') ? (char)(c - 'A' + 'a') : c; }  // 1 文字版
inline char to_upper(char c) { return ('a' <= c && c <= 'z') ? (char)(c - 'a' + 'A') : c; }
inline char swap_case(char c) { return to_lower(c) != c ? to_lower(c) : to_upper(c); }
template <class Seq>
str_internal::seq_t<Seq> caesar_shift(const Seq& s0, long long k, const Alphabet& A = Alphabet::lower()) {  // A の文字を番号で k ずらす（負も。A にない文字はそのまま）
    auto s = str_internal::seq_t<Seq>(str_internal::seq(s0));
    const long long sg = A.size(), sh = ((k % sg) + sg) % sg;
    for (auto& c : s)
        if (A.contains(c)) c = (std::decay_t<decltype(c)>)A.value((int)((A.id(c) + sh) % sg));
    return s;
}
inline char caesar_shift(char c, long long k, const Alphabet& A = Alphabet::lower()) {  // 1 文字版（A にない文字はそのまま）
    if (!A.contains(c)) return c;
    const long long sg = A.size(), sh = ((k % sg) + sg) % sg;
    return (char)A.value((int)((A.id(c) + sh) % sg));
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
template <class A, class B, class W>
int common_prefix_length(const A& a0, const B& b0, W wild) {  // wild（'?' など）はどちらにあっても一致とする（ABC287 D）
    const auto& a = str_internal::seq(a0);
    const auto& b = str_internal::seq(b0);
    const long long w = str_internal::code_of(wild);
    const std::size_t k = a.size() < b.size() ? a.size() : b.size();
    std::size_t i = 0;
    for (; i < k; ++i) {
        const long long x = str_internal::code_of(a[i]), y = str_internal::code_of(b[i]);
        if (x != y && x != w && y != w) break;
    }
    return (int)i;
}
template <class A, class B, class W>
int common_suffix_length(const A& a0, const B& b0, W wild) {
    const auto& a = str_internal::seq(a0);
    const auto& b = str_internal::seq(b0);
    const long long w = str_internal::code_of(wild);
    const std::size_t k = a.size() < b.size() ? a.size() : b.size();
    std::size_t i = 0;
    for (; i < k; ++i) {
        const long long x = str_internal::code_of(a[a.size() - 1 - i]), y = str_internal::code_of(b[b.size() - 1 - i]);
        if (x != y && x != w && y != w) break;
    }
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
template <class T, class S>
std::vector<int> subsequence_positions(const T& t0, const S& s0) {  // t を s の部分列として左から貪欲に当てはめた位置（部分列でなければ空。t が空でも空なので size() == |t| で判定）  O(|s|)
    const auto& t = str_internal::seq(t0);
    const auto& s = str_internal::seq(s0);
    std::vector<int> r;
    for (std::size_t i = 0; i < s.size() && r.size() < t.size(); ++i)
        if (str_internal::code_of(s[i]) == str_internal::code_of(t[r.size()])) r.push_back((int)i);
    if (r.size() < t.size()) r.clear();
    return r;
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
inline std::vector<std::string> split(const std::string& s, const std::string& sep, bool skip_empty = false) {  // 文字列 sep で区切る（", " など）
    if (sep.empty()) str_internal::fail("split", "区切りが空");
    std::vector<std::string> r;
    std::size_t b = 0;
    while (true) {
        const std::size_t i = s.find(sep, b);
        const std::size_t e = i == std::string::npos ? s.size() : i;
        if (!(skip_empty && e == b)) r.push_back(s.substr(b, e - b));
        if (i == std::string::npos) break;
        b = i + sep.size();
    }
    return r;
}
inline std::string strip_left(const std::string& s, const std::string& chars = " \t\r\n") {  // 先頭から chars の文字を除く（既定は空白・タブ・改行）
    const std::size_t i = s.find_first_not_of(chars);
    return i == std::string::npos ? std::string() : s.substr(i);
}
inline std::string strip_right(const std::string& s, const std::string& chars = " \t\r\n") {  // 末尾から chars の文字を除く（ABC366 B・ABC367 B）
    const std::size_t i = s.find_last_not_of(chars);
    return i == std::string::npos ? std::string() : s.substr(0, i + 1);
}
inline std::string strip(const std::string& s, const std::string& chars = " \t\r\n") { return strip_left(strip_right(s, chars), chars); }  // 両端から
inline std::string strip_left(const std::string& s, char c) { return strip_left(s, std::string(1, c)); }
inline std::string strip_right(const std::string& s, char c) { return strip_right(s, std::string(1, c)); }
inline std::string strip(const std::string& s, char c) { return strip(s, std::string(1, c)); }
template <class Cont>
std::string join(const Cont& vs, const std::string& sep = " ") {  // 区切りを挟んでつなげる（vector・set・deque など。数も文字列にする。小数は fixed で 15 桁。答えの出力に）
    using V = std::decay_t<decltype(*std::begin(vs))>;
    std::string r;
    bool first = true;
    for (const auto& x : vs) {
        if (!first) r += sep;
        first = false;
        if constexpr (std::is_same<V, std::string>::value) r += x;
        else if constexpr (std::is_same<V, char>::value) r += x;
        else if constexpr (std::is_constructible<std::string, const V&>::value) r += std::string(x);  // const char*・string_view
        else if constexpr (std::is_integral<V>::value) r += std::to_string(x);
        else {
            std::ostringstream os;
            if constexpr (std::is_floating_point<V>::value) os << std::fixed << std::setprecision(15);
            os << x;
            r += os.str();
        }
    }
    return r;
}
template <class Cont>
std::string join(const Cont& vs, char sep) { return join(vs, std::string(1, sep)); }  // 区切りを 1 文字で（'\n' など）

// ---- 隣どうしの条件で区切る・隣が同じにならない並べ替え ----
template <class Seq, class F>
std::vector<StrSpan> split_runs_by(const Seq& s0, F same) {  // same(s[i-1], s[i]) が成り立つ間を 1 つにまとめた区間（ABC456 C：隣が違う区間）  O(N)
    const auto& s = str_internal::seq(s0);
    std::vector<StrSpan> r;
    const int n = (int)s.size();
    for (int i = 0; i < n;) {
        int j = i + 1;
        while (j < n && same(s[j - 1], s[j])) ++j;
        r.push_back({i, j});
        i = j;
    }
    return r;
}
namespace str_internal {
template <class Seq>
seq_t<Seq> rearrange_impl(const Seq& s0, bool lex_smallest) {  // 隣が同じにならない並べ替え（できなければ空）
    const auto& s = seq(s0);
    using E = elem_t<Seq>;
    const long long n = (long long)s.size();
    std::map<long long, long long> cnt;  // 値 → 個数
    std::map<long long, E> rep;          // 値 → 元の要素
    for (std::size_t i = 0; i < s.size(); ++i) {
        const long long v = code_of(s[i]);
        ++cnt[v], rep.emplace(v, s[i]);
    }
    seq_t<Seq> res = seq_t<Seq>(s);
    long long mx = 0;
    for (auto& [v, c] : cnt) mx = std::max(mx, c);
    if (2 * mx > n + 1) return seq_t<Seq>();
    if (!lex_smallest) {  // 多い順に偶数番目 → 奇数番目に詰める
        std::vector<std::pair<long long, long long>> order;
        for (auto& [v, c] : cnt) order.push_back({-c, v});
        std::sort(order.begin(), order.end());
        long long pos = 0;
        for (auto& [nc, v] : order)
            for (long long k = 0; k < -nc; ++k) {
                res[(std::size_t)pos] = rep[v];
                pos += 2;
                if (pos >= n) pos = 1;
            }
        return res;
    }
    // 辞書順最小：残り rem が奇数で個数 (rem+1)/2 の値があればそれ、なければ直前と違う最小の値
    std::set<std::pair<long long, long long>> by_cnt;  // (個数, 値)
    for (auto& [v, c] : cnt) by_cnt.insert({c, v});
    bool has_prev = false;
    long long prev = 0;
    for (long long i = 0; i < n; ++i) {
        const long long rem = n - i;
        long long pick;
        const auto top = *by_cnt.rbegin();
        if (rem % 2 == 1 && top.first == (rem + 1) / 2) pick = top.second;
        else {
            auto it = cnt.begin();
            if (has_prev && it->first == prev) ++it;
            pick = it->first;
        }
        res[(std::size_t)i] = rep[pick];
        auto it = cnt.find(pick);
        by_cnt.erase({it->second, pick});
        if (--it->second == 0) cnt.erase(it);
        else by_cnt.insert({it->second, pick});
        has_prev = true, prev = pick;
    }
    return res;
}
}  // namespace str_internal
template <class Seq>
str_internal::seq_t<Seq> rearrange_no_adjacent_equal(const Seq& s) {  // 隣どうしが違うように並べ替えた 1 つ（できなければ空。ABC459 D）  O(N log σ)
    return str_internal::rearrange_impl(s, false);
}
template <class Seq>
str_internal::seq_t<Seq> lex_smallest_no_adjacent_equal(const Seq& s) {  // 隣どうしが違う並べ替えのうち辞書順最小（できなければ空）  O(N log σ)
    return str_internal::rearrange_impl(s, true);
}
#endif  // GG_STR_BASE
#ifndef GG_PALINDROME
#define GG_PALINDROME
// ---- 回文（文字列第二世代）  Manacher・区間の回文判定 O(1)・最長の回文・回文の個数・回文にする・回文部分列・回文分割・s 以下の回文の数 ----
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
template <class T, class Seq>
T count_palindromes_leq(const Seq& s0, const Alphabet& A = Alphabet::lower()) {  // A の文字からなる長さ |s| の回文で s 以下（Alphabet の順）のものの数（ABC242 E）  O(N)
    const auto& s = str_internal::seq(s0);
    const int n = (int)s.size(), h = (n + 1) / 2;
    const char* func = "count_palindromes_leq";
    T r = T(0);
    for (int i = 0; i < h; ++i) r = str_internal::checked_add(str_internal::checked_mul(r, T(A.size()), func), T(A.id(s[i])), func);  // 前半が小さいもの
    int cmp = 0;  // 前半をそのまま折り返した回文と s の比較
    for (int i = h; i < n && cmp == 0; ++i) {
        const int a = A.id(s[n - 1 - i]), b = A.id(s[i]);
        cmp = a < b ? -1 : (a > b ? 1 : 0);
    }
    if (cmp <= 0) r = str_internal::checked_add(r, T(1), func);
    return r;
}
#endif  // GG_PALINDROME
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト
// =============================================================================
namespace pal_test {

#define PAL_CHECK(...)                                                                          \
    do {                                                                                        \
        if (!(__VA_ARGS__)) {                                                                   \
            std::cerr << "[palindrome] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
            std::exit(1);                                                                       \
        }                                                                                       \
    } while (0)

#if defined(_DEBUG) || defined(_GLIBCXX_DEBUG) || (defined(__GNUC__) && !defined(__OPTIMIZE__))
[[maybe_unused]] static const int TEST_SCALE = 1;
#else
[[maybe_unused]] static const int TEST_SCALE = 8;
#endif

[[maybe_unused]] static long long rnd(std::mt19937_64& rng, long long lo, long long hi) {  // [lo, hi] の一様乱数（hi - lo が 2^64 未満）
    const unsigned long long w = (unsigned long long)hi - (unsigned long long)lo + 1;
    return w == 0 ? (long long)rng() : (long long)((unsigned long long)lo + rng() % w);
}
// 大きさがばらばらの long long（境目の値を多めに）
[[maybe_unused]] static long long rnd_ll(std::mt19937_64& rng) {
    static const long long special[] = {0, 1, -1, 2, -2, 3, -3, 7, -7, 10, -10, 2147483647LL, -2147483647LL, -2147483648LL, 2147483648LL,
                                        4294967295LL, 4294967296LL, -4294967296LL, 1000000000LL, -1000000000LL, 1000000007LL, 998244353LL,
                                        3037000499LL, 3037000500LL, -3037000500LL, 1000000000000000000LL, -1000000000000000000LL,
                                        4611686018427387904LL, -4611686018427387904LL, 9223372036854775807LL, -9223372036854775807LL - 1,
                                        9223372036854775806LL, -9223372036854775807LL};
    const int t = (int)(rng() % 8);
    if (t == 0) return special[rng() % (sizeof(special) / sizeof(special[0]))];
    const int bits = (int)(rng() % 64);
    long long v = (long long)(rng() >> (63 - bits) >> 1);
    if (t <= 1) v = (long long)(rng() % 21) - 10;
    return (rng() & 1) ? v : -v;
}

// ---- 文字列のテスト用（文字列第二世代） ----
// 'base' から sigma 種類の文字で長さ n のランダムな文字列（形をいろいろに：一様・同じ文字の連続・周期・フィボナッチ・Thue–Morse・ほぼ同じ文字）
[[maybe_unused]] static std::string rnd_str(std::mt19937_64& rng, int n, int sigma, char base = 'a') {
    std::string s(n, base);
    for (int i = 0; i < n; ++i) s[i] = (char)(base + (int)(rng() % (unsigned)sigma));
    return s;
}
[[maybe_unused]] static std::string rnd_str_shaped(std::mt19937_64& rng, int n, int sigma, char base = 'a') {
    std::string s(n, base);
    const int t = (int)(rng() % 7);
    if (t == 0 || sigma == 1) return rnd_str(rng, n, sigma, base);
    if (t == 1) {  // 同じ文字の連続
        for (int i = 0; i < n;) {
            const char c = (char)(base + (int)(rng() % (unsigned)sigma));
            int len = 1 + (int)(rng() % 6);
            while (len-- > 0 && i < n) s[i++] = c;
        }
    } else if (t == 2) {  // 周期
        const int p = 1 + (int)(rng() % 5);
        std::string u = rnd_str(rng, p, sigma, base);
        for (int i = 0; i < n; ++i) s[i] = u[i % p];
        if (n > 0 && rng() % 2) s[rng() % n] = (char)(base + (int)(rng() % (unsigned)sigma));
    } else if (t == 3) {  // フィボナッチ文字列
        std::string a = std::string(1, base), b = std::string(1, (char)(base + 1)) + base;
        while ((int)b.size() < n) {
            std::string c = b + a;
            a = b, b = c;
        }
        s = b.substr(0, n);
    } else if (t == 4) {  // Thue–Morse 列
        for (int i = 0; i < n; ++i) s[i] = (char)(base + (__builtin_popcount((unsigned)i) & 1));
    } else if (t == 5) {  // ほぼ同じ文字
        for (int i = 0; i < n; ++i) s[i] = (rng() % 8 == 0) ? (char)(base + (int)(rng() % (unsigned)sigma)) : base;
    } else {
        return rnd_str(rng, n, sigma, base);
    }
    return s;
}

static bool bf_pal(const std::string& s, int l, int r) {
    for (int i = l, j = r - 1; i < j; ++i, --j)
        if (s[i] != s[j]) return false;
    return true;
}

static void test_basic(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 3000 * TEST_SCALE; ++it) {
        const int n = (int)(rng() % 16), sg = 1 + (int)(rng() % 3);
        const std::string s = rnd_str_shaped(rng, n, sg);
        const Palindromes P(s);
        PAL_CHECK(P.size() == n);
        long long cnt = 0;
        StrSpan best{0, 0};
        for (int l = 0; l <= n; ++l)
            for (int r = l; r <= n; ++r) {
                const bool p = bf_pal(s, l, r);
                PAL_CHECK(P.is_palindrome(l, r) == p);
                if (p && r > l) {
                    ++cnt;
                    if (r - l > best.len()) best = {l, r};
                }
            }
        PAL_CHECK(P.count_all() == cnt && P.longest() == best);
        const auto L = center_palindrome_lengths(s);
        PAL_CHECK((int)L.size() == std::max(0, 2 * n - 1));
        for (int c = 0; c < 2 * n - 1; ++c) {
            int mx = 0;
            for (int l = 0; l <= n; ++l) {
                const int r = c + 1 - l;
                if (r >= l && r <= n && r - l >= mx && bf_pal(s, l, r) && (c % 2 == 0 ? r > l : true)) mx = std::max(mx, r - l);
            }
            PAL_CHECK(L[c] == mx && P.length_at_center(c) == mx && P.span_at_center(c).len() == mx);
        }
        for (int i = 0; i < n; ++i) PAL_CHECK(2 * P.radius_odd(i) - 1 == L[2 * i]);
        for (int i = 0; i <= n; ++i) PAL_CHECK(2 * P.radius_even(i) == (i == 0 || i == n ? 0 : L[2 * i - 1]));
        int pre = 0, suf = 0;
        for (int k = 1; k <= n; ++k) {
            if (bf_pal(s, 0, k)) pre = k;
            if (bf_pal(s, n - k, n)) suf = k;
        }
        PAL_CHECK(P.longest_prefix_palindrome() == pre && P.longest_suffix_palindrome() == suf);
        const auto st = P.longest_starting_at(), en = P.longest_ending_at();
        for (int i = 0; i < n; ++i) {
            int a = 0, b = 0;
            for (int k = 1; i + k <= n; ++k)
                if (bf_pal(s, i, i + k)) a = k;
            for (int k = 1; k <= i + 1; ++k)
                if (bf_pal(s, i + 1 - k, i + 1)) b = k;
            PAL_CHECK(st[i] == a && en[i] == b);
        }
        // 回文にする
        const std::string wp = min_palindrome_with_prefix(s), ws = min_palindrome_with_suffix(s);
        PAL_CHECK(is_palindrome(wp) && wp.compare(0, n, s) == 0 && (int)wp.size() == 2 * n - suf);
        PAL_CHECK(is_palindrome(ws) && ws.compare(ws.size() - n, n, s) == 0 && (int)ws.size() == 2 * n - pre);
        int ch = 0;
        for (int i = 0; i < n / 2; ++i) ch += s[i] != s[n - 1 - i];
        PAL_CHECK(min_changes_to_palindrome(s) == ch);
        std::map<char, int> mc;
        for (char x : s) ++mc[x];
        int odd = 0;
        for (auto& [k, v] : mc) odd += v & 1;
        PAL_CHECK(can_rearrange_to_palindrome(s) == (odd <= 1));
        // 回文分割（O(N^3) の DP）
        std::vector<int> dp(n + 1, 1 << 30);
        dp[0] = 0;
        for (int r = 1; r <= n; ++r)
            for (int l = 0; l < r; ++l)
                if (bf_pal(s, l, r)) dp[r] = std::min(dp[r], dp[l] + 1);
        const auto parts = palindrome_partition_min(s);
        PAL_CHECK((int)parts.size() == dp[n]);
        int pos = 0;
        for (auto sp : parts) {
            PAL_CHECK(sp.l == pos && sp.r > sp.l && bf_pal(s, sp.l, sp.r));
            pos = sp.r;
        }
        PAL_CHECK(pos == n);
        ++cases;
    }
}

static void test_subseq(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 2000 * TEST_SCALE; ++it) {
        const int n = (int)(rng() % 13), sg = 1 + (int)(rng() % 3);
        const std::string s = rnd_str_shaped(rng, n, sg);
        int best = 0;
        long long cnt = 0;
        for (int mask = 1; mask < (1 << n); ++mask) {
            std::string t;
            for (int i = 0; i < n; ++i)
                if (mask >> i & 1) t += s[i];
            if (is_palindrome(t)) ++cnt, best = std::max(best, (int)t.size());
        }
        PAL_CHECK(longest_palindromic_subsequence_length(s) == best);
        const std::string lps = longest_palindromic_subsequence(s);
        PAL_CHECK((int)lps.size() == best && is_palindrome(lps) && is_subsequence(lps, s));
        PAL_CHECK(count_palindromic_subsequences<long long>(s) == cnt && count_palindromic_subsequences<double>(s) == (double)cnt);
        // 数の列でも
        const std::vector<long long> v(s.begin(), s.end());
        PAL_CHECK(longest_palindromic_subsequence_length(v) == best && Palindromes(v).count_all() == Palindromes(s).count_all());
        ++cases;
    }
}


// 長さ n の全部の文字列（A の番号の列）を辞書順で
static std::vector<std::vector<int>> all_id_strings(int n, int sg) {
    std::vector<std::vector<int>> r(1);
    for (int k = 0; k < n; ++k) {
        std::vector<std::vector<int>> nr;
        for (const auto& v : r)
            for (int c = 0; c < sg; ++c) {
                auto w = v;
                w.push_back(c);
                nr.push_back(w);
            }
        r = nr;
    }
    return r;
}
static void test_final_additions(std::mt19937_64& rng, int& cases) {  // count_palindromes_leq
    for (int it = 0; it < 1500 * TEST_SCALE; ++it) {
        const int sg = 1 + (int)(rng() % 3), n = (int)(rng() % 7);
        std::string letters = std::string("abc").substr(0, sg);
        std::shuffle(letters.begin(), letters.end(), rng);
        const Alphabet A = Alphabet::of(letters);  // 独自の順も
        std::string s(n, ' ');
        for (char& c : s) c = letters[rng() % sg];
        const std::vector<int> sid = to_ids(s, A);
        long long cnt = 0;
        for (const auto& v : all_id_strings(n, sg)) cnt += std::equal(v.begin(), v.end(), v.rbegin()) && v <= sid;
        PAL_CHECK(count_palindromes_leq<long long>(s, A) == cnt && count_palindromes_leq<double>(s, A) == (double)cnt);
        ++cases;
    }
    PAL_CHECK(count_palindromes_leq<long long>(std::string(26, 'z')) == 2481152873203736576LL);  // 26^13
    PAL_CHECK(count_palindromes_leq<long long>("AXA", Alphabet::upper()) == 24 && count_palindromes_leq<long long>("ABCZAZ", Alphabet::upper()) == 29);  // ABC242 E
}
static void random_test() {
    std::mt19937_64 rng(20261012);
    int cases = 0;
    test_basic(rng, cases);
    test_subseq(rng, cases);
    test_final_additions(rng, cases);
    // 答えが long long に収まれば、途中の和が大きくても止まらない
    {
        const std::string s = "b" + std::string(15, 'a') + "b" + std::string(20, 'a') + "b" + std::string(7, 'a') + "b" + std::string(8, 'a') + "b" + std::string(12, 'a');
        PAL_CHECK(count_palindromic_subsequences<long long>(s) == 6172198434589173569LL);
        ++cases;
    }
    std::cout << "[palindrome] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[palindrome] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    auto ms = [](std::chrono::steady_clock::time_point a) {
        return (long long)std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count();
    };
    const int n = 10000000;
    const std::string s = rnd_str_shaped(rng, n, 2), a(n, 'a');
    auto t0 = std::chrono::steady_clock::now();
    const Palindromes P(s), PA(a);
    const long long t_build = ms(t0);
    t0 = std::chrono::steady_clock::now();
    long long x = 0;
    for (int i = 0; i < 10000000; ++i) {
        int l = (int)(rng() % (n + 1)), r = (int)(rng() % (n + 1));
        if (l > r) std::swap(l, r);
        x += P.is_palindrome(l, r);
    }
    const long long t_q = ms(t0);
    t0 = std::chrono::steady_clock::now();
    x += P.count_all() + PA.longest().len() + P.longest_starting_at()[3];
    const long long t_misc = ms(t0);
    const std::string m = rnd_str(rng, 5000, 3);
    t0 = std::chrono::steady_clock::now();
    x += (long long)longest_palindromic_subsequence(m).size() + (long long)palindrome_partition_min(m).size();
    const long long t_n2 = ms(t0);
    std::cout << "[palindrome] benchmark: Palindromes 構築（N = 10^7 × 2）" << t_build << " ms, is_palindrome × 10^7 " << t_q
              << " ms, count_all + longest + longest_starting_at（N = 10^7）" << t_misc << " ms, 最長回文部分列 + 回文分割（N = 5000）" << t_n2
              << " ms  (check " << x % 7 << ")\n";
}

}  // namespace pal_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace pal_demo {

[[maybe_unused]] static std::string to_s(bool b) { return b ? "true" : "false"; }
[[maybe_unused]] static std::string to_s(long double v) {
    if (std::fabs(v) < 5e-13L) v = 0;
    std::ostringstream os;
    os << std::setprecision(12) << v;
    return os.str();
}
[[maybe_unused]] static std::string to_s(double v) { return to_s((long double)v); }
[[maybe_unused]] static std::string to_s(const std::string& s) { return s; }
[[maybe_unused]] static std::string to_s(const char* s) { return s; }
template <class T>
[[maybe_unused]] static std::string to_s(const T& v) {  // << で出せるもの（整数・Int128・modint など）
    std::ostringstream os;
    os << v;
    return os.str();
}
template <class T>
[[maybe_unused]] static std::string to_s(const std::vector<T>& v) {
    std::string s = "[";
    for (size_t i = 0; i < v.size(); ++i) s += (i ? ", " : "") + to_s(v[i]);
    return s + "]";
}
[[maybe_unused]] static std::string pad(const std::string& s, size_t w) {
    size_t cols = 0;
    for (size_t i = 0; i < s.size(); ++i) {
        const unsigned char c = (unsigned char)s[i];
        if ((c & 0xC0) == 0x80) continue;
        cols += (c >= 0xE0 ? 2 : 1);
    }
    return cols >= w ? s : s + std::string(w - cols, ' ');
}
[[maybe_unused]] static void banner(const std::string& s) {
    std::cout << "\n=====================================================================\n"
              << "  " << s << "\n"
              << "=====================================================================\n";
}
[[maybe_unused]] static void api(const std::string& usage, const std::string& ret, const std::string& what, const std::string& cost) {
    std::cout << "\n--------------------------------------------------------------------\n"
              << "■ " << usage << "\n"
              << "    返り値 : " << ret << "\n"
              << "    内容   : " << what << "\n"
              << "    計算量 : " << cost << "\n";
}
[[maybe_unused]] static void head(const std::string& s) {
    std::cout << "\n--------------------------------------------------------------------\n"
              << "■ " << s << "\n";
}
[[maybe_unused]] static void note(const std::string& s) { std::cout << "    " << s << "\n"; }
[[maybe_unused]] static void code(const std::string& s) { std::cout << "  " << s << "\n"; }
[[maybe_unused]] static void stmt(const std::string& s) { code(s + ";"); }
template <class V>
[[maybe_unused]] static void show(const char* expr, const V& v, const std::string& comment) {
    std::string line = "  " + pad(expr, 34) + " -> " + to_s(v);
    if (!comment.empty()) line = pad(line, 58) + "  // " + comment;
    std::cout << line << "\n";
}
#define SHOW(expr, comment) pal_demo::show(#expr, (expr), comment)
#define RUN(...)                   \
    pal_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("回文（palindrome） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  Palindromes P(s)：is_palindrome(l, r) / radius_odd / radius_even / length_at_center / span_at_center / longest / count_all\n"
              << "                    longest_prefix_palindrome / longest_suffix_palindrome / longest_starting_at / longest_ending_at\n"
              << "  center_palindrome_lengths(s)                         Manacher の生の値（中心 2N-1 個）\n"
              << "  min_palindrome_with_prefix / min_palindrome_with_suffix   s で始まる / 終わる最短の回文\n"
              << "  min_changes_to_palindrome / can_rearrange_to_palindrome  書き換え・並べ替えて回文に\n"
              << "  longest_palindromic_subsequence(_length) / palindrome_partition_min / count_palindromic_subsequences<T>\n"
              << "  count_palindromes_leq<T>(s, A)                       同じ長さで s 以下の回文の数\n";

    api("【1】 Palindromes P(s) / P.is_palindrome(l, r)", "bool", "s[l, r) が回文か（Manacher で前処理）。何度でも O(1)", "構築 O(N)・判定 O(1)");
    RUN(const std::string S = "abacabad");
    RUN(const Palindromes P(S));
    SHOW(P.is_palindrome(0, 7), "\"abacaba\"");
    SHOW(P.is_palindrome(1, 4), "\"bac\"");
    note("問題の「L 文字目から R 文字目」は P.is_palindrome(L - 1, R)。全体だけなら is_palindrome(s)（00_str_base）。");

    api("【2】 P.radius_odd(i) / P.radius_even(i) / P.length_at_center(c) / center_palindrome_lengths(s)", "int / vector<int>",
        "s[i] が中心の最長の回文の半径 r（長さ 2r - 1）/ s[i-1] と s[i] の間が中心の半径 r（長さ 2r）/ 中心 c = 0〜2N-2 ごとの長さ", "O(1) / O(N)");
    SHOW(P.radius_odd(3), "\"abacaba\" は半径 4");
    SHOW(P.radius_even(3), "\"a|c\" の間は回文なし");
    SHOW(center_palindrome_lengths("abba"), "c = 0〜6（奇数の c が文字の間）");

    api("【3】 P.longest() / P.count_all() / P.longest_prefix_palindrome() / P.longest_suffix_palindrome()", "StrSpan / long long / int",
        "最長の回文（同じ長さなら左）・回文である部分文字列の数（位置が違えば別）・最長の回文の接頭辞 / 接尾辞", "O(N)");
    SHOW(substr_of(S, P.longest()), "");
    SHOW(P.count_all(), "");
    SHOW(P.longest_prefix_palindrome(), "\"abacaba\"");
    SHOW(P.longest_suffix_palindrome(), "\"d\"");

    api("【4】 P.longest_starting_at() / P.longest_ending_at()", "vector<int>",
        "各位置から始まる / で終わる最長の回文の長さ（回文で区切る DP の前処理などに）", "O(N)");
    SHOW(P.longest_starting_at(), "");
    SHOW(P.longest_ending_at(), "");

    api("【5】 min_palindrome_with_prefix(s) / min_palindrome_with_suffix(s)", "string",
        "s で始まる最短の回文（後ろに足す）/ s で終わる最短の回文（前に足す）", "O(N)");
    SHOW(min_palindrome_with_prefix("ABC"), "ABC398 F（ABCBA）の入力例 1");
    SHOW(min_palindrome_with_prefix("abcba"), "もう回文");
    SHOW(min_palindrome_with_suffix("aab"), "前に足す");

    api("【6】 min_changes_to_palindrome(s) / can_rearrange_to_palindrome(s)", "int / bool",
        "1 文字ずつ書き換えて回文にする最小の回数 / 並べ替えて回文にできるか（奇数個の文字が 1 種類以下）", "O(N)");
    SHOW(min_changes_to_palindrome("abcd"), "a ↔ d、b ↔ c");
    SHOW(can_rearrange_to_palindrome("aabbc"), "abcba");
    SHOW(can_rearrange_to_palindrome("aabc"), "");

    api("【7】 longest_palindromic_subsequence(s) / _length(s) / count_palindromic_subsequences<T>(s)", "string / int / T",
        "最長回文部分列（飛び飛び）を 1 つ・その長さ・回文である部分列の数（位置が違えば別。T は mint や long long）", "O(N^2)");
    SHOW(longest_palindromic_subsequence("character"), "");
    SHOW(longest_palindromic_subsequence_length("character"), "");
    SHOW(count_palindromic_subsequences<long long>("aab"), "a, a, b, aa");
    note("答えを mod で出すときは math_v2 の using mint = ModInt998; で count_palindromic_subsequences<mint>(s)。");
    note("long long で溢れるときはメッセージを出して止まる。");

    api("【8】 palindrome_partition_min(s)", "vector<StrSpan>", "回文に分ける最小の個数の分け方（.size() が個数）", "O(N^2)");
    RUN(const auto parts = palindrome_partition_min("abaccd"));
    SHOW(parts.size(), "aba | cc | d");
    SHOW(substr_of("abaccd", parts[1]), "");
    note("N が大きい（10^5 以上）ときや、異なる回文の種類・出現回数は 04_palindromic_tree。");

    api("【9】 count_palindromes_leq<T>(s, A)", "T（mint・long long。long long で溢れたら止める）",
        "A の文字からなる長さ |s| の回文のうち、辞書順（A の番号の順）で s 以下のものの数（前半を σ 進数で読む）", "O(N)");
    SHOW(count_palindromes_leq<long long>("AXA", Alphabet::upper()), "ABC242 E の入力例 1（AAA〜AXA の 24 個）");
    SHOW(count_palindromes_leq<long long>("ABCZAZ", Alphabet::upper()), "");
    note("答えを mod で出すときは count_palindromes_leq<mint>(S, Alphabet::upper())。同じ長さの文字列全体での順位は 05_lex_order の lex_rank_same_length<T>。");
    note("k 文字まで書き換えてよい回文（中心ごと）は 01_rolling_hash の center_palindrome_lengths_with_changes（ABC468 D）。");
}

}  // namespace pal_demo

// =============================================================================
//  Windows のコンソールで日本語が文字化けしないようにする（Windows 以外では何もしない）
// =============================================================================
#ifdef _WIN32
extern "C" __declspec(dllimport) int __stdcall SetConsoleOutputCP(unsigned int);
#endif
static void setup_console() {
#ifdef _WIN32
    const char* probe = "あ";  // UTF-8 なら 0xE3 0x81 0x82、CP932 なら 0x82 0xA0
    SetConsoleOutputCP((unsigned char)probe[0] == 0xE3 ? 65001u : 932u);
#endif
}

int main() {
    setup_console();
    pal_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    pal_test::benchmark();    // 速度計測
    pal_demo::demo();         // 使い方デモ
    return 0;
}
