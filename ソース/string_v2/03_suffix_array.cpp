// =============================================================================
//  接尾辞配列（SA-IS・LCP・O(1) の LCP と比較・出現・異なる部分文字列・k 番目・LCP の和・複数の文字列・最長共通部分文字列）   --- 文字列第二世代
// -----------------------------------------------------------------------------
//  分類       : 03 接尾辞
//  依存       : 00_str_base（同梱済み）
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）
//  文字コード : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  列は string でも vector<int>・vector<long long> でもよい（数は中で圧縮する）。位置は 0 始まり、部分文字列は [l, r)
//  ACL の suffix_array / lcp_array とは名前を変えている（make_suffix_array / make_lcp_array）
//  「部分文字列」は連続する部分。飛び飛びの「部分列」の LCS は 06_lcs_edit
//
//  機能一覧                                                                       計算量
//    make_suffix_array(s) / make_lcp_array(s, sa)   接尾辞配列（SA-IS）/ 隣どうしの LCP（Kasai。長さ N - 1）   O(N)
//    SuffixArray SA(s)                     sa / rank / lcp_adj（隣どうしの LCP）を持つ                  O(N)
//    SuffixArray::of_strings(vs)           複数の文字列をまとめて（区切りの位置も sa に入る。owner = -1）
//      SA.owner(p) / offset(p) / position(t, off) / length_at(p)   何番目の文字列か / その中の位置 / 逆 / 接尾辞の長さ
//      SA.lcp(i, j)                        位置 i, j から始まる接尾辞の LCP                              O(1)（初回 O(N log N)）
//      SA.compare(l1, r1, l2, r2)          部分文字列の辞書順（-1 / 0 / 1）                               O(1)
//      SA.find_range(pat) / count(pat) / find_all(pat)   pat で始まる SA の区間 / 出現回数 / 出現位置     O(M log N)
//      SA.count_occurrences(l, r) / occurrence_range(l, r)   s[l, r) の出現回数 / SA の区間            O(log N)
//      SA.count_distinct_substrings()      異なる部分文字列（空でない）の数                               O(N)
//      SA.distinct_substrings_by_length()  長さごとの異なる部分文字列の数                                 O(N)
//      SA.kth_substring(k)                 異なる部分文字列の辞書順 k 番目（0 から。StrSpan）              O(log N)
//      SA.substring_rank(l, r)             異なる部分文字列のうち s[l, r) より小さいものの数               O(log N)
//      SA.kth_substring_with_dup(k) / count_substrings_with_dup()   位置が違えば別として k 番目 / 総数     O(log N)
//      SA.longest_repeated()               2 か所以上に現れる最長（重なってよい。CommonSubstring{len, pos_a, pos_b}）  O(N)
//      SA.longest_repeated_non_overlapping()   重ならずに 2 か所に現れる最長                              O(N log N)
//      SA.sum_lcp_with_all() / sum_lcp_all_pairs()   各 i で Σ_j lcp(i, j)（j = i も）/ Σ_{i<j} lcp   O(N)
//    longest_common_substring(a, b)        2 つの最長共通部分文字列（CommonSubstring{len, pos_a, pos_b}）  O(N + M)
//    longest_common_substrings_all(a, b)   最長共通部分文字列をすべて（a での位置。辞書順・重複なし）      O(N + M)
//    longest_common_substring_k(vs)        K 個すべてに現れる最長の部分文字列（vs[0] での StrSpan）        O(全体)
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    sa_test::random_test()  … 全部の部分文字列を並べた愚直解と全機能を突き合わせるランダムテスト
//    sa_test::benchmark()    … 速度計測
//    sa_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
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
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト
// =============================================================================
namespace sa_test {

#define SA_CHECK(...)                                                                          \
    do {                                                                                        \
        if (!(__VA_ARGS__)) {                                                                   \
            std::cerr << "[suffix_array] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
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

using VL = std::vector<long long>;
// 愚直：文字列の組（of_strings と同じ並び：各文字列のあとに区切り）
struct Brute {
    std::vector<VL> strs;
    VL all;             // 値（区切りは使わない）
    std::vector<int> owner, lim;
    std::vector<VL> keys;  // 位置 p の接尾辞の並べ替え用の列（区切りは文字より小さく互いに違う）
    int n = 0;
    bool is_multi = false;
    explicit Brute(const std::vector<VL>& vs, bool multi) : strs(vs), is_multi(multi) {
        const int K = (int)vs.size();
        std::vector<long long> code;  // 並べ替え用の値
        for (int t = 0; t < K; ++t) {
            const int st = (int)all.size();
            for (long long x : vs[t]) all.push_back(x), owner.push_back(t), lim.push_back(st + (int)vs[t].size()), code.push_back(x);
            if (multi) all.push_back(0), owner.push_back(-1), lim.push_back((int)all.size() - 1), code.push_back(std::numeric_limits<long long>::min() + t);
        }
        n = (int)all.size();
        for (int p = 0; p < n; ++p) keys.push_back(VL(code.begin() + p, code.end()));
    }
    std::vector<int> sa() const {
        std::vector<int> r(n);
        std::iota(r.begin(), r.end(), 0);
        std::sort(r.begin(), r.end(), [&](int a, int b) { return keys[a] < keys[b]; });
        return r;
    }
    VL sub(int l, int r) const { return VL(all.begin() + l, all.begin() + r); }
    int lcp(int i, int j) const {
        if (i == n || j == n) return 0;
        int k = 0;
        while (i + k < lim[i] && j + k < lim[j] && all[i + k] == all[j + k]) ++k;
        return k;
    }
};

static void check_sa(std::mt19937_64& rng, const SuffixArray& S, const Brute& B, int& cases) {
    const int n = B.n;
    SA_CHECK(S.size() == n && S.sa == B.sa());
    for (int k = 0; k < n; ++k) SA_CHECK(S.rank[S.sa[k]] == k && S.length_at(k) == B.lim[k] - k);
    for (int k = 0; k + 1 < n; ++k) SA_CHECK(S.lcp_adj[k] == B.lcp(S.sa[k], S.sa[k + 1]));
    for (int i = 0; i <= n; ++i)
        for (int j = 0; j <= n; ++j) SA_CHECK(S.lcp(i, j) == B.lcp(i, j));
    // 部分文字列の一覧
    std::vector<std::pair<int, int>> subs;  // (l, r)
    for (int l = 0; l < n; ++l)
        for (int r = l + 1; r <= B.lim[l]; ++r) subs.push_back({l, r});
    std::vector<VL> multi;
    for (auto [l, r] : subs) multi.push_back(B.sub(l, r));
    std::sort(multi.begin(), multi.end());
    std::vector<VL> dist = multi;
    dist.erase(std::unique(dist.begin(), dist.end()), dist.end());
    SA_CHECK(S.count_distinct_substrings() == (long long)dist.size() && S.count_substrings_with_dup() == (long long)multi.size());
    for (long long k = 0; k < (long long)dist.size(); ++k) {
        const StrSpan sp = S.kth_substring(k);
        SA_CHECK(B.sub(sp.l, sp.r) == dist[k]);
    }
    for (long long k = 0; k < (long long)multi.size(); ++k) {
        const StrSpan sp = S.kth_substring_with_dup(k);
        SA_CHECK(B.sub(sp.l, sp.r) == multi[k]);
    }
    const auto byl = S.distinct_substrings_by_length();
    std::vector<long long> exl(byl.size(), 0);
    for (const auto& v : dist) {
        SA_CHECK(v.size() < byl.size());
        ++exl[v.size()];
    }
    SA_CHECK(byl == exl);
    for (auto [l, r] : subs) {
        const VL v = B.sub(l, r);
        long long occ = 0;
        for (auto [l2, r2] : subs)
            if (r2 - l2 == r - l && B.sub(l2, r2) == v) ++occ;
        SA_CHECK(S.count_occurrences(l, r) == occ);
        SA_CHECK(S.substring_rank(l, r) == (long long)(std::lower_bound(dist.begin(), dist.end(), v) - dist.begin()));
        ++cases;
    }
    for (int l = 0; l <= n; ++l) SA_CHECK(S.count_occurrences(l, l) == n + (B.is_multi ? 0 : 1));  // 空の部分文字列
    // 比較
    for (int t = 0; t < 60 && !subs.empty(); ++t) {
        auto [l1, r1] = subs[rng() % subs.size()];
        auto [l2, r2] = subs[rng() % subs.size()];
        if (rng() % 4 == 0) r1 = l1;
        const VL a = B.sub(l1, r1), b = B.sub(l2, r2);
        SA_CHECK(S.compare(l1, r1, l2, r2) == (a < b ? -1 : (a == b ? 0 : 1)));
    }
    // パターン
    for (int t = 0; t < 30; ++t) {
        VL pat;
        if (!subs.empty() && rng() % 2) {
            auto [l, r] = subs[rng() % subs.size()];
            pat = B.sub(l, r);
            if (rng() % 3 == 0 && !pat.empty() && pat.back() > -1000000 && pat.back() < 1000000) pat.back() += (rng() % 2 ? 1 : -1);
        } else {
            for (int k = (int)(rng() % 4); k > 0; --k) pat.push_back(B.all.empty() ? 0 : B.all[rng() % B.all.size()]);
        }
        std::vector<int> ex;
        for (int p = 0; p < n; ++p)
            if (p + (int)pat.size() <= B.lim[p] && B.sub(p, p + (int)pat.size()) == pat) ex.push_back(p);
        if (pat.empty()) {  // 空なら置ける位置すべて（1 つの文字列なら 0〜N、of_strings なら 0〜N-1 = 各文字列の 0〜長さ）
            ex.clear();
            for (int p = 0; p < n + (B.is_multi ? 0 : 1); ++p) ex.push_back(p);
        }
        SA_CHECK(S.find_all(pat) == ex && S.count(pat) == (long long)ex.size());
    }
    // 繰り返し
    int best = 0, best_nov = 0;
    for (auto [l, r] : subs)
        for (int q = l + 1; q < n; ++q)
            if (q + (r - l) <= B.lim[q] && B.sub(q, q + (r - l)) == B.sub(l, r)) {
                best = std::max(best, r - l);
                if (q >= r) best_nov = std::max(best_nov, r - l);
            }
    const CommonSubstring c = S.longest_repeated(), d = S.longest_repeated_non_overlapping();
    SA_CHECK(c.len == best && d.len == best_nov);
    if (best) {
        SA_CHECK(c.pos_a < c.pos_b && B.sub(c.pos_a, c.pos_a + best) == B.sub(c.pos_b, c.pos_b + best) && c.pos_b + best <= B.lim[c.pos_b]);
        VL smallest;  // 最長のうち辞書順最小
        bool any = false;
        for (auto [l, r] : subs)
            if (r - l == best && S.count_occurrences(l, r) >= 2 && (!any || B.sub(l, r) < smallest)) smallest = B.sub(l, r), any = true;
        SA_CHECK(B.sub(c.pos_a, c.pos_a + best) == smallest);
    } else {
        SA_CHECK(c.pos_a == -1 && c.pos_b == -1);
    }
    if (best_nov) SA_CHECK(d.pos_a + best_nov <= d.pos_b && B.sub(d.pos_a, d.pos_a + best_nov) == B.sub(d.pos_b, d.pos_b + best_nov));
    // LCP の和
    const auto sw = S.sum_lcp_with_all();
    long long pairs = 0;
    for (int i = 0; i < n; ++i) {
        long long t = 0;
        for (int j = 0; j < n; ++j) t += i == j ? B.lim[i] - i : B.lcp(i, j);
        for (int j = i + 1; j < n; ++j) pairs += B.lcp(i, j);
        SA_CHECK(sw[i] == t);
    }
    SA_CHECK(S.sum_lcp_all_pairs() == pairs);
}

static void test_single(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 500 * TEST_SCALE; ++it) {
        const int n = (int)(rng() % 13), sg = 1 + (int)(rng() % 3);
        const std::string s = rnd_str_shaped(rng, n, sg);
        const SuffixArray S(s);
        check_sa(rng, S, Brute({VL(s.begin(), s.end())}, false), cases);
        SA_CHECK(make_suffix_array(s) == S.sa && make_lcp_array(s, S.sa) == S.lcp_adj);
        // 数の列（大きな値・負の値）
        VL a(n);
        const long long vals[] = {std::numeric_limits<long long>::min(), -5, 0, 1000000000000000000LL, std::numeric_limits<long long>::max()};
        for (auto& x : a) x = vals[rng() % (2 + sg)];
        check_sa(rng, SuffixArray(a), Brute({a}, false), cases);
    }
}

static void test_multi(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 300 * TEST_SCALE; ++it) {
        const int K = 1 + (int)(rng() % 3), sg = 1 + (int)(rng() % 3);
        std::vector<std::string> vs(K);
        std::vector<VL> vv(K);
        for (int t = 0; t < K; ++t) vs[t] = rnd_str_shaped(rng, (int)(rng() % 6), sg), vv[t] = VL(vs[t].begin(), vs[t].end());
        const SuffixArray S = SuffixArray::of_strings(vs);
        const Brute B(vv, true);
        check_sa(rng, S, B, cases);
        for (int p = 0; p < B.n; ++p) {
            SA_CHECK(S.owner(p) == B.owner[p]);
            if (B.owner[p] >= 0) SA_CHECK(S.position(S.owner(p), S.offset(p)) == p);
        }
        // 2 つ・K 個の最長共通部分文字列
        if (K >= 2) {
            const std::string& a = vs[0];
            const std::string& b = vs[1];
            int best = 0;
            std::set<std::string> alls;
            for (int i = 0; i < (int)a.size(); ++i)
                for (int j = 0; j < (int)b.size(); ++j) {
                    int k = 0;
                    while (i + k < (int)a.size() && j + k < (int)b.size() && a[i + k] == b[j + k]) ++k;
                    if (k > best) best = k, alls.clear();
                    if (k == best && k > 0) alls.insert(a.substr(i, k));
                }
            const CommonSubstring c = longest_common_substring(a, b);
            SA_CHECK(c.len == best);
            if (best) SA_CHECK(a.substr(c.pos_a, best) == b.substr(c.pos_b, best));
            else SA_CHECK(c.pos_a == -1);
            std::vector<std::string> got;
            for (StrSpan sp : longest_common_substrings_all(a, b)) got.push_back(substr_of(a, sp));
            SA_CHECK(got == std::vector<std::string>(alls.begin(), alls.end()));
            // K 個
            int bk = 0;
            for (int i = 0; i < (int)a.size(); ++i)
                for (int L = 1; i + L <= (int)a.size(); ++L) {
                    bool ok = true;
                    for (int t = 1; t < K; ++t) ok &= vs[t].find(a.substr(i, L)) != std::string::npos;
                    if (ok) bk = std::max(bk, L);
                }
            const StrSpan sk = longest_common_substring_k(vs);
            SA_CHECK(sk.len() == bk);
            for (int t = 1; t < K && bk; ++t) SA_CHECK(vs[t].find(substr_of(a, sk)) != std::string::npos);
        }
    }
}

static void test_sais(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 40 * TEST_SCALE; ++it) {
        const int n = (int)(rng() % (it % 8 == 0 ? 3000 : 400)), sg = 1 + (int)(rng() % 4);
        const std::string s = rnd_str_shaped(rng, n, sg);
        std::vector<int> ex(n);
        std::iota(ex.begin(), ex.end(), 0);
        std::sort(ex.begin(), ex.end(), [&](int a, int b) { return s.compare(a, std::string::npos, s, b, std::string::npos) < 0; });
        const auto sa = make_suffix_array(s);
        SA_CHECK(sa == ex);
        const auto lcp = make_lcp_array(s, sa);
        for (int k = 0; k + 1 < n; ++k) {
            int h = 0;
            while (sa[k] + h < n && sa[k + 1] + h < n && s[sa[k] + h] == s[sa[k + 1] + h]) ++h;
            SA_CHECK(lcp[k] == h);
        }
        ++cases;
    }
}

static void random_test() {
    std::mt19937_64 rng(20261012);
    int cases = 0;
    test_single(rng, cases);
    test_multi(rng, cases);
    test_sais(rng, cases);
    std::cout << "[suffix_array] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[suffix_array] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    auto ms = [](std::chrono::steady_clock::time_point a) {
        return (long long)std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count();
    };
    const int n = 1000000;
    const std::string r26 = rnd_str(rng, n, 26), r2 = rnd_str_shaped(rng, n, 2), aa(n, 'a');
    auto t0 = std::chrono::steady_clock::now();
    const auto s1 = make_suffix_array(r26);
    const auto s2 = make_suffix_array(r2);
    const auto s3 = make_suffix_array(aa);
    const long long t_sa = ms(t0);
    t0 = std::chrono::steady_clock::now();
    const auto l1 = make_lcp_array(r26, s1);
    const long long t_lcp = ms(t0);
    t0 = std::chrono::steady_clock::now();
    const SuffixArray S(r2);
    const long long t_build = ms(t0);
    t0 = std::chrono::steady_clock::now();
    long long x = 0;
    for (int i = 0; i < 1000000; ++i) x += S.lcp((int)(rng() % n), (int)(rng() % n));
    const long long t_q = ms(t0);
    t0 = std::chrono::steady_clock::now();
    for (int i = 0; i < 100000; ++i) x += S.kth_substring_with_dup((long long)(rng() % (unsigned long long)S.count_substrings_with_dup())).len();
    const long long t_k = ms(t0);
    std::vector<long long> v(n);
    for (auto& y : v) y = rnd_ll(rng);
    t0 = std::chrono::steady_clock::now();
    const auto s4 = make_suffix_array(v);
    const long long t_v = ms(t0);
    std::cout << "[suffix_array] benchmark: make_suffix_array（N = 10^6。σ = 26 / 2 / 全部 a）合計 " << t_sa << " ms, make_lcp_array " << t_lcp
              << " ms, SuffixArray 構築 " << t_build << " ms, lcp(i, j) × 10^6（表の構築込み）" << t_q << " ms, kth_substring_with_dup × 10^5（初回の準備込み）"
              << t_k << " ms, long long の列（N = 10^6）" << t_v << " ms  (check " << (s1[7] + s2[5] + s3[3] + l1[9] + x + s4[1]) % 7 << ")\n";
}

}  // namespace sa_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace sa_demo {

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
#define SHOW(expr, comment) sa_demo::show(#expr, (expr), comment)
#define RUN(...)                   \
    sa_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("接尾辞配列（suffix_array） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  make_suffix_array / make_lcp_array                   接尾辞配列と隣どうしの LCP（ACL の suffix_array / lcp_array と同じ）\n"
              << "  SuffixArray SA(s)：sa / rank / lcp_adj / lcp(i, j) / compare / find_range / count / find_all / count_occurrences\n"
              << "     count_distinct_substrings / distinct_substrings_by_length / kth_substring / substring_rank / kth_substring_with_dup\n"
              << "     longest_repeated / longest_repeated_non_overlapping / sum_lcp_with_all / sum_lcp_all_pairs\n"
              << "  SuffixArray::of_strings(vs)：複数の文字列（owner / offset / position）\n"
              << "  longest_common_substring / longest_common_substrings_all / longest_common_substring_k   最長共通部分文字列（連続）\n";

    api("【1】 make_suffix_array(s) / make_lcp_array(s, sa)", "vector<int>",
        "接尾辞を辞書順に並べた開始位置（SA-IS）と、隣どうしの最長共通接頭辞（長さ N - 1）。数の列にも", "O(N)");
    RUN(const std::string S = "banana");
    SHOW(make_suffix_array(S), "a, ana, anana, banana, na, nana");
    SHOW(make_lcp_array(S, make_suffix_array(S)), "");

    api("【2】 SuffixArray SA(s) / SA.lcp(i, j) / SA.compare(l1, r1, l2, r2)", "int",
        "任意の 2 つの接尾辞の LCP・部分文字列の辞書順を O(1) で（初回に疎テーブルを作る）", "O(1)");
    RUN(const SuffixArray SA(S));
    SHOW(SA.sa, "");
    SHOW(SA.rank, "rank[p] = 位置 p の順位");
    SHOW(SA.lcp(1, 3), "\"anana\" と \"ana\"");
    SHOW(SA.compare(0, 2, 2, 4), "\"ba\" と \"na\"");

    api("【3】 SA.count(pat) / SA.find_all(pat) / SA.count_occurrences(l, r)", "long long / vector<int>",
        "パターンの出現回数・位置（何度でも O(M log N)）/ 部分文字列 s[l, r) の出現回数（O(log N)）", "—");
    SHOW(SA.count("ana"), "重なりも");
    SHOW(SA.find_all("an"), "");
    SHOW(SA.count_occurrences(1, 3), "\"an\"");
    note("ABC362 G（Count Substring Query）: クエリごとに SA.count(T_i)。");

    api("【4】 count_distinct_substrings / distinct_substrings_by_length / kth_substring(k) / substring_rank(l, r)", "long long / vector / StrSpan",
        "異なる部分文字列の数・長さごとの数・辞書順 k 番目（0 から）・s[l, r) が何番目か", "O(N) / O(log N)");
    SHOW(SA.count_distinct_substrings(), "a, an, ana, … の 15 個");
    SHOW(SA.distinct_substrings_by_length(), "r[len]");
    SHOW(substr_of(S, SA.kth_substring(0)), "辞書順で最初");
    SHOW(substr_of(S, SA.kth_substring(4)), "a, an, ana, anan, anana");
    SHOW(SA.substring_rank(1, 4), "\"ana\" より小さい異なる部分文字列の数");

    api("【5】 kth_substring_with_dup(k) / count_substrings_with_dup()", "StrSpan / long long",
        "位置が違えば別として数えたときの辞書順 k 番目（同じ文字列が出現回数だけ並ぶ）。複数の文字列でも（ABC280 Ex）", "O(log N)（初回 O(N)）");
    SHOW(SA.count_substrings_with_dup(), "N(N+1)/2");
    SHOW(substr_of(S, SA.kth_substring_with_dup(2)), "a, a, a, an, …");
    SHOW(substr_of(S, SA.kth_substring_with_dup(3)), "");

    api("【6】 longest_repeated() / longest_repeated_non_overlapping()", "CommonSubstring{len, pos_a, pos_b}",
        "2 か所以上に現れる最長の部分文字列（重なってよい / 重ならない）。ないときは len = 0", "O(N) / O(N log N)");
    SHOW(SA.longest_repeated().len, "\"ana\"（位置 1 と 3。重なる）");
    SHOW(SA.longest_repeated_non_overlapping().len, "\"an\" か \"na\"");
    note("ABC141 E（Who Says a Pun?）は longest_repeated_non_overlapping().len。");

    api("【7】 sum_lcp_with_all() / sum_lcp_all_pairs()", "vector<long long> / long long",
        "各 i で Σ_j lcp(i, j)（j = i も含む）/ すべての組の LCP の和", "O(N)");
    SHOW(SA.sum_lcp_with_all(), "ABC213 F（Common Prefixes）の答えそのもの");
    SHOW(SA.sum_lcp_all_pairs(), "");

    api("【8】 SuffixArray::of_strings(vs) / owner(p) / offset(p) / position(t, off)", "SuffixArray",
        "複数の文字列をまとめて（区切りはどの文字より小さい。部分文字列は 1 つの文字列の中だけ）。sa には区切りの位置（owner = -1）も入る", "O(全体)");
    RUN(const SuffixArray M = SuffixArray::of_strings(std::vector<std::string>{"ab", "b"}));
    SHOW(M.sa, "位置 2 と 4 は区切り");
    SHOW(M.owner(3), "位置 3 は 1 番目の文字列");
    SHOW(M.count("b"), "両方の文字列での出現");
    SHOW(M.count_distinct_substrings(), "a, ab, b（文字列をまたがない）");

    api("【9】 longest_common_substring(a, b) / longest_common_substrings_all(a, b) / longest_common_substring_k(vs)", "CommonSubstring / vector<StrSpan> / StrSpan",
        "最長共通部分文字列（連続する部分）。2 つ・すべて・K 個。飛び飛びの「部分列」の LCS は 06_lcs_edit", "O(全体)");
    RUN(const CommonSubstring c = longest_common_substring("xabcde", "abcyabc"));
    SHOW(c.len, "");
    SHOW(c.pos_a, "a での位置");
    SHOW(c.pos_b, "b での位置");
    SHOW(longest_common_substrings_all("abxcd", "cdyab").size(), "\"ab\" と \"cd\"");
    SHOW(substr_of("abcd", longest_common_substring_k(std::vector<std::string>{"abcd", "zbcx", "bcq"})), "3 つすべてに");
}

}  // namespace sa_demo

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
    sa_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    sa_test::benchmark();    // 速度計測
    sa_demo::demo();         // 使い方デモ
    return 0;
}
