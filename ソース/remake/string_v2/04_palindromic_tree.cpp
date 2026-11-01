// =============================================================================
//  回文木（eertree：異なる回文・出現回数・各接頭辞の最長の回文接尾辞・回文への最小分割と分割数（N log N）・回文の包含）   --- 文字列第二世代
// -----------------------------------------------------------------------------
//  分類       : 04 回文
//  依存       : 00_str_base（同梱済み）
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）
//  文字コード : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  異なる回文は N 個以下で、それぞれが 1 つのノード。ノード 0 は長さ -1 の根、ノード 1 は空文字列、2 以降が回文
//  add(c) で 1 文字ずつ足せる（オンライン）。遷移は ノード数 × σ 個の int。区間の回文判定・最長の回文は 04_palindrome
//  回文の包含：回文 v の中の回文は link(v)（最長の回文接尾辞）と parent(v)（両端を除いたもの）をたどると全部出てくる（ABC237 Ex）
//
//  機能一覧                                                                       計算量
//    PalindromicTree pt(s, A) / pt.add(c)  作る / 1 文字足す（最長の回文接尾辞のノードを返す）              O(N σ) / O(1)（ならし）
//    pt.count_distinct() / count_total()   異なる回文の数 / 回文の部分文字列（位置が違えば別）の数           O(1) / O(N)
//    pt.occurrences()                      ノードごとの出現回数                                           O(N)
//    pt.palindromes() / span(v)            異なる回文の最初の出現（StrSpan）                               O(N) / O(1)
//    pt.min_partition() / min_partition_spans()   回文に分ける最小の個数 / 分け方                          O(N log N)
//    pt.count_partitions<T>()              回文に分ける方法の数（T は mint・long long）                     O(N log N)
//    ノード: size() / length() / len(v) / link(v) / parent(v) / next(v, c) / node_of_prefix(i) / first_end(v)
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    pt_test::random_test()  … 全部の部分文字列を調べる愚直解・O(N^2) の DP と全機能を突き合わせるランダムテスト
//    pt_test::benchmark()    … 速度計測
//    pt_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
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
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト
// =============================================================================
namespace pt_test {

#define PT_CHECK(...)                                                                          \
    do {                                                                                        \
        if (!(__VA_ARGS__)) {                                                                   \
            std::cerr << "[palindromic_tree] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
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

struct Mod {  // テスト用の mod 998244353
    unsigned v = 0;
    Mod() = default;
    Mod(long long x) : v((unsigned)(((x % 998244353) + 998244353) % 998244353)) {}
    Mod operator+(const Mod& o) const { return Mod((long long)v + o.v); }
    bool operator==(const Mod& o) const { return v == o.v; }
};
static bool is_pal(const std::string& s) { return std::equal(s.begin(), s.end(), s.rbegin()); }

static void random_test() {
    std::mt19937_64 rng(20261012);
    int cases = 0;
    for (int it = 0; it < 2000 * TEST_SCALE; ++it) {
        const int sg = 1 + (int)(rng() % 3);
        const std::string s = rnd_str_shaped(rng, (int)(rng() % 30), sg);
        const int n = (int)s.size();
        const PalindromicTree pt(s);
        // 異なる回文と出現回数
        std::map<std::string, long long> occ;
        std::map<std::string, int> first_end;
        for (int r = 1; r <= n; ++r)
            for (int l = 0; l < r; ++l) {
                const std::string t = s.substr(l, r - l);
                if (!is_pal(t)) continue;
                ++occ[t];
                if (!first_end.count(t)) first_end[t] = r - 1;
            }
        PT_CHECK(pt.count_distinct() == (int)occ.size() && pt.size() == (int)occ.size() + 2 && pt.length() == n);
        long long total = 0;
        for (auto& [t, c] : occ) total += c;
        PT_CHECK(pt.count_total() == total);
        const auto oc = pt.occurrences();
        const auto pals = pt.palindromes();
        std::set<std::string> seen;
        for (int v = 2; v < pt.size(); ++v) {
            const StrSpan sp = pt.span(v);
            const std::string t = s.substr(sp.l, sp.len());
            PT_CHECK(pals[v - 2] == sp && is_pal(t) && (int)t.size() == pt.len(v) && occ.count(t) && !seen.count(t));
            seen.insert(t);
            PT_CHECK(oc[v] == occ[t] && pt.first_end(v) == first_end[t]);
            // link：最長の真の回文接尾辞 / parent：両端を除いたもの / next：両側に足したもの
            int best = 0;
            for (int k = 1; k < (int)t.size(); ++k)
                if (is_pal(t.substr(t.size() - k))) best = k;
            PT_CHECK(pt.len(pt.link(v)) == best);
            const int pv = pt.parent(v);
            PT_CHECK((int)t.size() <= 2 ? pv == 1 : (pt.len(pv) == (int)t.size() - 2 && s.substr(pt.span(pv).l, pt.len(pv)) == t.substr(1, t.size() - 2)));
            if ((int)t.size() >= 3) PT_CHECK(pt.next(pv, t[0]) == v);
        }
        // 各接頭辞の最長の回文接尾辞
        for (int i = 0; i <= n; ++i) {
            int best = 0;
            for (int k = 1; k <= i; ++k)
                if (is_pal(s.substr(i - k, k))) best = k;
            PT_CHECK(pt.len(pt.node_of_prefix(i)) == best);
        }
        // 分割
        std::vector<int> dp(n + 1, 1 << 29);
        std::vector<long long> ways(n + 1, 0);
        dp[0] = 0, ways[0] = 1;
        for (int i = 1; i <= n; ++i)
            for (int j = 0; j < i; ++j)
                if (is_pal(s.substr(j, i - j))) dp[i] = std::min(dp[i], dp[j] + 1), ways[i] += ways[j];
        PT_CHECK(pt.min_partition() == dp[n] && pt.count_partitions<long long>() == ways[n] && pt.count_partitions<Mod>() == Mod(ways[n]));
        const auto mpp = pt.min_partition_prefix();
        const auto cnp = pt.count_partitions_prefix<long long>();
        PT_CHECK((int)mpp.size() == n + 1 && (int)cnp.size() == n + 1);
        for (int i = 0; i <= n; ++i) PT_CHECK(mpp[i] == dp[i] && cnp[i] == ways[i]);
        const auto parts = pt.min_partition_spans();
        int pos = 0;
        for (const StrSpan& sp : parts) {
            PT_CHECK(sp.l == pos && is_pal(s.substr(sp.l, sp.len())));
            pos = sp.r;
        }
        PT_CHECK(pos == n && (int)parts.size() == dp[n]);
        ++cases;
    }
    // 1 文字ずつ足しても同じ・数の列
    for (int it = 0; it < 300 * TEST_SCALE; ++it) {
        std::vector<int> a(rng() % 25);
        for (auto& x : a) x = (int)(rng() % 3);
        PalindromicTree pt(Alphabet::ints(3));
        for (int x : a) pt.add(x);
        std::set<std::vector<int>> st;
        for (int l = 0; l < (int)a.size(); ++l)
            for (int r = l + 1; r <= (int)a.size(); ++r) {
                std::vector<int> t(a.begin() + l, a.begin() + r);
                if (std::equal(t.begin(), t.end(), t.rbegin())) st.insert(t);
            }
        PT_CHECK(pt.count_distinct() == (int)st.size());
        ++cases;
    }
    // 問い合わせと add を混ぜる（出現回数の作り直し）
    for (int it = 0; it < 300 * TEST_SCALE; ++it) {
        const int sg = 1 + (int)(rng() % 3);
        PalindromicTree pt;
        std::string s;
        for (int k = 0; k < 15; ++k) {
            const char c = (char)('a' + rng() % sg);
            s += c, pt.add(c);
            const auto& oc = pt.occurrences();
            long long total = 0;
            for (int v = 2; v < pt.size(); ++v) {
                const std::string t = s.substr(pt.span(v).l, pt.len(v));
                long long c2 = 0;
                for (std::size_t l = 0; l + t.size() <= s.size(); ++l) c2 += s.compare(l, t.size(), t) == 0;
                PT_CHECK(oc[v] == c2);
                total += c2;
            }
            PT_CHECK(pt.count_total() == total && pt.min_partition_prefix().back() == pt.min_partition());
        }
        ++cases;
    }
    std::cout << "[palindromic_tree] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[palindromic_tree] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    auto ms = [](std::chrono::steady_clock::time_point a) {
        return (long long)std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count();
    };
    const std::string s = rnd_str(rng, 1000000, 2);
    auto t0 = std::chrono::steady_clock::now();
    const PalindromicTree pt(s, Alphabet::of("ab"));
    const long long t_build = ms(t0);
    t0 = std::chrono::steady_clock::now();
    long long x = pt.count_total() + pt.count_distinct();
    const long long t_cnt = ms(t0);
    t0 = std::chrono::steady_clock::now();
    x += pt.min_partition();
    const long long t_part = ms(t0);
    const std::string a = std::string(1000000, 'a');
    t0 = std::chrono::steady_clock::now();
    const PalindromicTree pa(a);
    x += pa.count_partitions<double>() > 0;
    const long long t_aa = ms(t0);
    std::cout << "[palindromic_tree] benchmark: 作成（N = 10^6, σ = 2）" << t_build << " ms, count_total " << t_cnt << " ms, min_partition " << t_part
              << " ms, 全部 a（10^6）の作成 + count_partitions " << t_aa << " ms  (check " << x % 7 << ")\n";
}

}  // namespace pt_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace pt_demo {

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
#define SHOW(expr, comment) pt_demo::show(#expr, (expr), comment)
#define RUN(...)                   \
    pt_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("回文木（palindromic_tree） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  add / count_distinct / count_total / occurrences / palindromes / span\n"
              << "  min_partition / min_partition_spans / min_partition_prefix / count_partitions<T> / count_partitions_prefix<T>\n"
              << "  size / length / len / link / parent / next / node_of_prefix / first_end   ノード（0 は長さ -1 の根、1 は空文字列）\n";

    RUN(const PalindromicTree pt(std::string("ababb")));

    api("【1】 count_distinct() / count_total() / occurrences() / palindromes()", "int / long long / vector<long long> / vector<StrSpan>",
        "異なる回文の数 / 回文の部分文字列の数（位置が違えば別）/ ノードごとの出現回数 / それぞれの最初の出現", "O(N)");
    SHOW(pt.count_distinct(), "a, b, aba, bab, bb");
    SHOW(pt.count_total(), "");
    SHOW(pt.palindromes(), "");
    SHOW(pt.occurrences(), "ノード 0・1 は 0");

    api("【2】 ノード：len(v) / link(v) / parent(v) / next(v, c) / node_of_prefix(i)", "int",
        "長さ / 最長の真の回文接尾辞 / 両端を除いた回文 / 両側に c を足した回文（なければ -1）/ s[0, i) の最長の回文接尾辞", "O(1)");
    {
        RUN(const int v = pt.node_of_prefix(4));
        SHOW(pt.len(v), "\"abab\" の最長の回文接尾辞は bab");
        SHOW(pt.len(pt.link(v)), "bab の最長の真の回文接尾辞は b");
        SHOW(pt.len(pt.parent(v)), "bab の両端を除くと a");
        SHOW(pt.next(pt.parent(v), 'b') == v, "");
    }

    api("【3】 min_partition() / min_partition_spans() / count_partitions<T>() / min_partition_prefix() / count_partitions_prefix<T>()",
        "int / vector<StrSpan> / T / vector<int> / vector<T>",
        "回文に分ける最小の個数 / その分け方 / 分け方の数（N = 10^6 でも速い。N ≤ 5000 なら 04_palindrome の palindrome_partition_min でもよい）", "O(N log N)");
    SHOW(pt.min_partition(), "aba + bb");
    SHOW(pt.min_partition_spans(), "");
    SHOW(pt.count_partitions<long long>(), "");
    SHOW(pt.min_partition_prefix(), "各接頭辞 s[0, i) の最小の個数（count_partitions_prefix<T>() は各接頭辞の分け方の数）");

    head("問題での使い方：ABC237 Ex（Hakata）");
    note("回文どうしが含まれない最大の選び方 → 回文のノードを頂点に、v → link(v)・v → parent(v) をたどれる関係（推移閉包）で");
    note("「どれも含まない」最大の集合（反鎖）= 頂点の数 − 二部グラフの最大マッチング（graph_grid_v2 08_bipartite_matching）。入力例:");
    for (const std::string s : {"ababb", "xyz", "xxxxxxxxxx"}) {
        const PalindromicTree t(s);
        const int m = t.size();
        std::vector<std::vector<char>> reach(m, std::vector<char>(m, 0));
        for (int v = 2; v < m; ++v) {  // v から link・parent をたどって着く回文（ノード 2 以上）
            std::vector<int> st = {v};
            while (!st.empty()) {
                const int u = st.back();
                st.pop_back();
                for (int w : {t.link(u), t.parent(u)})
                    if (w >= 2 && !reach[v][w]) reach[v][w] = 1, st.push_back(w);
            }
        }
        // 二部マッチング（左 v → 右 w：v が w を含む）
        std::vector<int> mt(m, -1);
        int match = 0;
        for (int v = 2; v < m; ++v) {
            std::vector<char> used(m, 0);
            auto dfs = [&](auto&& self, int x) -> bool {
                for (int w = 2; w < m; ++w)
                    if (reach[x][w] && !used[w]) {
                        used[w] = 1;
                        if (mt[w] < 0 || self(self, mt[w])) return mt[w] = x, true;
                    }
                return false;
            };
            match += dfs(dfs, v);
        }
        note("  " + s + " → " + to_s(t.count_distinct() - match));
    }
}

}  // namespace pt_demo

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
    pt_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    pt_test::benchmark();    // 速度計測
    pt_demo::demo();         // 使い方デモ
    return 0;
}
