// =============================================================================
//  接尾辞オートマトン（部分文字列か・出現回数・最初の出現・異なる部分文字列の数・辞書順 k 番目・最長共通部分文字列・状態（DP 用））   --- 文字列第二世代
// -----------------------------------------------------------------------------
//  分類       : 03 接尾辞
//  依存       : 00_str_base（同梱済み）
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）
//  文字コード : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  SuffixAutomaton sam(s, A)：s のすべての部分文字列を受け付ける最小の DFA。状態は 2N 個以下（0 が根 = 空文字列）
//  遷移は 状態数 × σ 個の int（N = 5×10^5・σ = 26 で約 100MB）。extend(c) で 1 文字ずつ足せる。add_string(t) で複数の文字列をまとめる
//  位置は足した文字の通し番号（1 つの文字列なら s の位置そのもの）。多くの問題は 03_suffix_array でも解ける（どちらか書きやすい方）
//
//  機能一覧                                                                       計算量
//    SuffixAutomaton sam(s, A) / sam.extend(c) / sam.add_string(t)   作る / 1 文字足す / 別の文字列として足す   O(N σ)
//    sam.contains(t) / count_occurrences(t) / first_occurrence(t)   部分文字列か / 出現回数 / 最初の開始位置   O(|t|)
//    sam.count_distinct_substrings() / total_length_distinct()     異なる部分文字列の数 / 長さの合計          O(状態数)
//    sam.kth_substring(k, with_dup) / kth_substring_ids(k, with_dup)   辞書順 k 番目（0 から。with_dup で重複も） O(|答え| σ)
//    sam.count_substrings_with_dup()       位置が違えば別として数えた数                                     O(状態数)
//    sam.longest_common_substring(t)       t の部分文字列でこれにも現れる最長（t での StrSpan）              O(|t|)
//    状態: size() / len(v) / link(v) / next(v, c) / next_id(v, id) / endpos_count(v) / first_end(v) / order_by_len() / last()
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    sam_test::random_test()  … 全部の部分文字列を数える愚直解と全機能を突き合わせるランダムテスト
//    sam_test::benchmark()    … 速度計測
//    sam_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
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
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト
// =============================================================================
namespace sam_test {

#define SAM_CHECK(...)                                                                          \
    do {                                                                                        \
        if (!(__VA_ARGS__)) {                                                                   \
            std::cerr << "[suffix_automaton] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
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

// 愚直：文字列の組（つながらない）の全部の部分文字列
struct Brute {
    std::vector<std::string> ss;
    std::map<std::string, long long> occ;      // 出現回数
    std::map<std::string, long long> first;    // 最初の開始位置（通し番号）
    explicit Brute(const std::vector<std::string>& v) : ss(v) {
        long long base = 0;
        for (const auto& s : ss) {
            for (int l = 0; l < (int)s.size(); ++l)
                for (int r = l + 1; r <= (int)s.size(); ++r) {
                    const std::string t = s.substr(l, r - l);
                    ++occ[t];
                    if (!first.count(t)) first[t] = base + l;
                    else first[t] = std::min(first[t], base + l);
                }
            base += (long long)s.size();
        }
    }
};

static void check(std::mt19937_64& rng, const SuffixAutomaton& sam, const Brute& B, int sg) {
    long long total_len = 0, dup = 0, chars = 0;
    for (const auto& s : B.ss) chars += (long long)s.size();
    for (auto& [t, c] : B.occ) total_len += (long long)t.size(), dup += c;
    SAM_CHECK(sam.count_distinct_substrings() == (long long)B.occ.size() && sam.total_length_distinct() == total_len && sam.count_substrings_with_dup() == dup);
    // 部分文字列か・回数・最初
    for (int it = 0; it < 30; ++it) {
        std::string t = rnd_str(rng, (int)(rng() % 6), sg + 1);
        if (!B.occ.empty() && rng() % 2) {
            auto itr = B.occ.begin();
            std::advance(itr, rng() % B.occ.size());
            t = itr->first;
        }
        const bool in = B.occ.count(t) > 0;
        SAM_CHECK(sam.contains(t) == (in || t.empty()));
        SAM_CHECK(sam.count_occurrences(t) == (t.empty() ? chars + (long long)std::max<std::size_t>(B.ss.size(), 1) : (in ? B.occ.at(t) : 0)));
        SAM_CHECK(sam.first_occurrence(t) == (t.empty() ? 0 : (in ? B.first.at(t) : -1)));
        // 最長共通部分文字列（t の中で）
        int best = 0, bl = 0;
        for (int l = 0; l < (int)t.size(); ++l)
            for (int r = l + 1; r <= (int)t.size(); ++r)
                if (B.occ.count(t.substr(l, r - l)) && r - l > best) best = r - l, bl = l;
        const StrSpan sp = sam.longest_common_substring(t);
        SAM_CHECK(sp.len() == best && (best == 0 || sp.l == bl));
    }
    // k 番目
    std::vector<std::string> dist, dups;
    for (auto& [t, c] : B.occ) {
        dist.push_back(t);
        for (long long k = 0; k < c; ++k) dups.push_back(t);
    }
    for (long long k = 0; k < (long long)dist.size(); ++k) SAM_CHECK(sam.kth_substring(k) == dist[k]);
    for (long long k = 0; k < (long long)dups.size(); ++k) SAM_CHECK(sam.kth_substring(k, true) == dups[k]);
    // 状態
    const auto ord = sam.order_by_len();
    SAM_CHECK((int)ord.size() == sam.size() && ord[0] == 0 && sam.link(0) == -1);
    for (int k = 1; k < sam.size(); ++k) SAM_CHECK(sam.len(ord[k - 1]) <= sam.len(ord[k]));
    for (int v = 1; v < sam.size(); ++v) SAM_CHECK(sam.len(sam.link(v)) < sam.len(v));
    // 状態ごとの回数：その状態の最長の文字列の回数
    for (auto& [t, c] : B.occ) {
        int v = 0;
        for (char ch : t) v = sam.next(v, ch);
        SAM_CHECK(v > 0 && sam.endpos_count(v) == c && sam.first_end(v) == B.first.at(t) + (long long)t.size() - 1);
        SAM_CHECK(sam.len(sam.link(v)) < (int)t.size() && (int)t.size() <= sam.len(v));
    }
}

static void random_test() {
    std::mt19937_64 rng(20261012);
    int cases = 0;
    for (int it = 0; it < 1500 * TEST_SCALE; ++it) {
        const int sg = 1 + (int)(rng() % 3);
        const std::string s = rnd_str_shaped(rng, (int)(rng() % 25), sg);
        const SuffixAutomaton sam(s);
        SAM_CHECK(s.size() < 2 || sam.size() <= 2 * (int)s.size() - 1);
        check(rng, sam, Brute({s}), sg);
        // 1 文字ずつ足しても同じ
        SuffixAutomaton s2;
        for (char c : s) s2.extend(c);
        SAM_CHECK(s2.size() == sam.size() && s2.count_distinct_substrings() == sam.count_distinct_substrings());
        ++cases;
    }
    // 複数の文字列
    for (int it = 0; it < 800 * TEST_SCALE; ++it) {
        const int sg = 1 + (int)(rng() % 3);
        std::vector<std::string> v(1 + rng() % 4);
        for (auto& s : v) s = rnd_str_shaped(rng, (int)(rng() % 10), sg);
        SuffixAutomaton sam;
        for (const auto& s : v) sam.add_string(s);
        check(rng, sam, Brute(v), sg);
        ++cases;
    }
    // extend・add_string と問い合わせを混ぜる（数え直し・last・続きの文字列）
    for (int it = 0; it < 300 * TEST_SCALE; ++it) {
        const int sg = 1 + (int)(rng() % 3);
        SuffixAutomaton sam;
        std::vector<std::string> v;
        for (int op = 0; op < 15; ++op) {
            if (v.empty() || rng() % 4 == 0) {
                v.push_back(rnd_str(rng, (int)(rng() % 4), sg));
                sam.add_string(v.back());
            } else {
                const char c = (char)('a' + rng() % sg);
                v.back() += c;
                sam.extend(c);
            }
            const Brute B(v);
            long long total_len = 0;
            for (auto& [t, c] : B.occ) total_len += (long long)t.size();
            SAM_CHECK(sam.count_distinct_substrings() == (long long)B.occ.size() && sam.total_length_distinct() == total_len);
            int w = 0;
            for (char ch : v.back()) w = sam.next(w, ch);
            SAM_CHECK(sam.last() == w);
            const std::string t = rnd_str(rng, (int)(rng() % 3), sg);
            long long chars = 0;
            for (auto& s : v) chars += (long long)s.size();
            SAM_CHECK(sam.count_occurrences(t) == (t.empty() ? chars + (long long)v.size() : (B.occ.count(t) ? B.occ.at(t) : 0)));
        }
        ++cases;
    }
    // 数の列
    for (int it = 0; it < 300 * TEST_SCALE; ++it) {
        std::vector<int> a(rng() % 20);
        for (auto& x : a) x = (int)(rng() % 3);
        const SuffixAutomaton sam(a, Alphabet::ints(3));
        std::set<std::vector<int>> st;
        for (int l = 0; l < (int)a.size(); ++l)
            for (int r = l + 1; r <= (int)a.size(); ++r) st.insert(std::vector<int>(a.begin() + l, a.begin() + r));
        SAM_CHECK(sam.count_distinct_substrings() == (long long)st.size());
        long long k = 0;
        for (const auto& t : st) SAM_CHECK(sam.kth_substring_ids(k++) == t);
        ++cases;
    }
    std::cout << "[suffix_automaton] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[suffix_automaton] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    auto ms = [](std::chrono::steady_clock::time_point a) {
        return (long long)std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count();
    };
    const std::string s = rnd_str(rng, 500000, 26);
    auto t0 = std::chrono::steady_clock::now();
    const SuffixAutomaton sam(s);
    const long long t_build = ms(t0);
    t0 = std::chrono::steady_clock::now();
    long long x = sam.count_distinct_substrings();
    for (int q = 0; q < 100000; ++q) {
        const int l = (int)(rng() % 499990);
        x += sam.count_occurrences(s.substr(l, 1 + rng() % 5));
    }
    const long long t_q = ms(t0);
    t0 = std::chrono::steady_clock::now();
    for (int q = 0; q < 100; ++q) x += (long long)sam.kth_substring((long long)(rng() % (unsigned long long)sam.count_substrings_with_dup()), true).size();
    const long long t_k = ms(t0);
    const std::string ab = rnd_str(rng, 1000000, 2);
    t0 = std::chrono::steady_clock::now();
    const SuffixAutomaton s2(ab, Alphabet::of("ab"));
    x += s2.size();
    const long long t_ab = ms(t0);
    std::cout << "[suffix_automaton] benchmark: 作成（N = 5×10^5, σ = 26）" << t_build << " ms, count_occurrences × 10^5（初回の数え込み）" << t_q
              << " ms, kth_substring(重複あり。答えの長さは平均 N / 3) × 100 " << t_k << " ms, 作成（N = 10^6, σ = 2）" << t_ab << " ms  (check " << x % 7 << ")\n";
}

}  // namespace sam_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace sam_demo {

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
#define SHOW(expr, comment) sam_demo::show(#expr, (expr), comment)
#define RUN(...)                   \
    sam_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("接尾辞オートマトン（suffix_automaton） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  contains / count_occurrences / first_occurrence / count_distinct_substrings / total_length_distinct\n"
              << "  kth_substring / kth_substring_ids / count_substrings_with_dup / longest_common_substring\n"
              << "  extend / add_string / size / len / link / next / next_id / endpos_count / first_end / order_by_len / last   状態（DP 用）\n";

    RUN(const SuffixAutomaton sam(std::string("missisippi")));

    api("【1】 contains(t) / count_occurrences(t) / first_occurrence(t)", "bool / long long / long long", "部分文字列か / 出現回数（重なりも）/ 最初の開始位置（なければ -1）", "O(|t|)");
    SHOW(sam.count_occurrences("i"), "ABC362 G の入力例 1：i が 4 回");
    SHOW(sam.count_occurrences("is"), "");
    SHOW(sam.contains("sip"), "");
    SHOW(sam.first_occurrence("ss"), "");

    api("【2】 count_distinct_substrings() / total_length_distinct() / count_substrings_with_dup()", "long long",
        "異なる部分文字列の数 / その長さの合計（この 2 つは足すたびに更新するので 1 文字ごとに聞いてもよい）/ 位置が違えば別として数えた数", "O(1) / O(1) / O(状態数 × σ)");
    SHOW(sam.count_distinct_substrings(), "");
    SHOW(sam.total_length_distinct(), "");
    SHOW(sam.count_substrings_with_dup(), "10 × 11 / 2");

    api("【3】 kth_substring(k, with_dup)", "string（kth_substring_ids は番号の列）", "辞書順 k 番目（0 から）。with_dup = true なら同じ文字列も出現の数だけ数える", "O(|答え| σ)");
    SHOW(sam.kth_substring(0), "");
    SHOW(sam.kth_substring(5), "i, ip, ipp, ippi, is, isi, ...");
    SHOW(sam.kth_substring(5, true), "i が 4 回、ip、ipp …");

    api("【4】 longest_common_substring(t)", "StrSpan（t での位置）", "t の部分文字列で、こちらにも現れる最長のもの（同じ長さなら左）", "O(|t|)");
    SHOW(sam.longest_common_substring("mississippi"), "\"missis\"");

    api("【5】 extend(c) / add_string(t)", "なし", "1 文字ずつ足す（最後に足した文字列の続き。オンライン）/ 別の文字列として足す（文字列の組の部分文字列をまとめて数える）", "O(σ)（ならし）");
    {
        RUN(SuffixAutomaton g);
        RUN(g.add_string("abab"), g.add_string("cab"));
        SHOW(g.count_distinct_substrings(), "abab と cab の部分文字列の種類");
        SHOW(g.count_occurrences("ab"), "2 つの文字列で合わせて 3 回");
        SHOW(g.first_occurrence("ca"), "通し番号（abab の 4 文字の後の 0 文字目）");
        SHOW(g.kth_substring(13, true), "ABC280 Ex の入力例 1 の 14 番目（2 番目の文字列の 1 文字目。どの文字列かは first_occurrence の通し番号から）");
        note("kth_substring は答えの長さに比例するので、答えが長くなる問い合わせを多数（ABC280 Ex の本番）なら 03_suffix_array の kth_substring_with_dup（StrSpan で返す。O(log N)）");
    }

    api("【6】 状態：len(v) / link(v) / next(v, c) / endpos_count(v) / order_by_len()", "int / int / int / long long / vector<int>",
        "状態 v は「出現位置の集合が同じ部分文字列」のまとまり。最長の長さ len、最短は len(link) + 1。order_by_len の逆順にたどると部分木の和（DP 用）",
        "O(1)（order_by_len は O(状態数 + N)）");
    {
        int v = 0;
        for (char c : std::string("ss")) v = sam.next(v, c);
        SHOW(sam.len(v), "\"ss\" の状態の最長（\"miss\"）");
        SHOW(sam.len(sam.link(v)), "最短は len(link) + 1 = 2");
        SHOW(sam.endpos_count(v), "");
        SHOW(sam.size(), "状態の数（2N 以下）");
    }
}

}  // namespace sam_demo

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
    sam_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    sam_test::benchmark();    // 速度計測
    sam_demo::demo();         // 使い方デモ
    return 0;
}
