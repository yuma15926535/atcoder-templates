// =============================================================================
//  最長増加部分列（LIS。長さ・復元・各位置で終わる / 始まる長さ・入りうる位置・数・減少・分割の最小数・2 次元）   --- 文字列第二世代
// -----------------------------------------------------------------------------
//  分類       : 06 部分列
//  依存       : 00_str_base（同梱済み）
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）
//  文字コード : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  列は vector<整数>（int・long long。負・10^18 級もよい）でも string でもよい
//  strict = true（既定）で狭義増加（a[i] < a[j]）、false で広義（a[i] ≤ a[j]。「減らない」）。添字は 0 始まり
//  数える関数は型 T を受け取る（math_v2 の using mint = ModInt998; で count_lis<mint>(a)。long long なら溢れたら止める）
//
//  機能一覧                                                                       計算量
//    lis_length(a, strict)                 最長増加部分列の長さ                                         O(N log N)
//    lis_restore(a, strict)                その 1 つ（添字の列。昇順）                                   O(N log N)
//    lis_ending_at(a, strict)              r[i] = a[i] で終わる最長の長さ                                 O(N log N)
//    lis_starting_at(a, strict)            r[i] = a[i] から始まる最長の長さ                               O(N log N)
//    lis_members(a, strict)                LisMembers{length, in_some, in_all}（入りうる・必ず入る位置）    O(N log N)
//    count_lis<T>(a, strict)               LIS の数（添字の列として）                                    O(N log N)
//    lds_length(a, strict)                 最長減少部分列の長さ                                          O(N log N)
//    min_increasing_cover(a, strict)       増加部分列に分ける最小の数（= 最長の減らない / 減る部分列）       O(N log N)
//    min_decreasing_cover(a, strict)       減少部分列に分ける最小の数                                    O(N log N)
//    lis_2d(ps)                            (x, y) が両方とも狭義に増える最長の列（添字。x の小さい順）       O(N log N)
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    lis_test::random_test()  … O(N^2) の DP・全部の部分列の列挙と全機能を突き合わせるランダムテスト
//    lis_test::benchmark()    … 速度計測
//    lis_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
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
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト
// =============================================================================
namespace lis_test {

#define LIS_CHECK(...)                                                                          \
    do {                                                                                        \
        if (!(__VA_ARGS__)) {                                                                   \
            std::cerr << "[lis] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
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

// テスト用の mod 998244353（数える関数に整数以外の型を渡す確認）
struct Mod {
    unsigned v = 0;
    Mod() = default;
    Mod(long long x) : v((unsigned)(((x % 998244353) + 998244353) % 998244353)) {}
    Mod operator+(const Mod& o) const { return Mod((long long)v + o.v); }
    Mod operator-(const Mod& o) const { return Mod((long long)v - o.v); }
    bool operator==(const Mod& o) const { return v == o.v; }
};

static bool less_by(long long x, long long y, bool strict) { return strict ? x < y : x <= y; }

// O(N^2) の DP：e[i] = i で終わる最長、s[i] = i から始まる最長、ce / cs = その数（mod）
struct Bf {
    std::vector<int> e, s;
    std::vector<Mod> ce, cs;
    int L = 0;
    Mod total;
};
static Bf bf_dp(const std::vector<long long>& v, bool strict) {
    const int n = (int)v.size();
    Bf b;
    b.e.assign(n, 1), b.s.assign(n, 1), b.ce.assign(n, Mod(1)), b.cs.assign(n, Mod(1));
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < i; ++j)
            if (less_by(v[j], v[i], strict)) {
                if (b.e[j] + 1 > b.e[i]) b.e[i] = b.e[j] + 1, b.ce[i] = b.ce[j];
                else if (b.e[j] + 1 == b.e[i]) b.ce[i] = b.ce[i] + b.ce[j];
            }
    for (int i = n - 1; i >= 0; --i)
        for (int j = i + 1; j < n; ++j)
            if (less_by(v[i], v[j], strict)) {
                if (b.s[j] + 1 > b.s[i]) b.s[i] = b.s[j] + 1, b.cs[i] = b.cs[j];
                else if (b.s[j] + 1 == b.s[i]) b.cs[i] = b.cs[i] + b.cs[j];
            }
    for (int i = 0; i < n; ++i) b.L = std::max(b.L, b.e[i]);
    b.total = n == 0 ? Mod(1) : Mod(0);
    for (int i = 0; i < n; ++i)
        if (b.e[i] == b.L) b.total = b.total + b.ce[i];
    return b;
}
// 増加部分列に分ける最小の数（バックトラック。n ≤ 9）
static int bf_cover(const std::vector<long long>& v, bool strict, bool increasing) {
    const int n = (int)v.size();
    for (int k = 0; k <= n; ++k) {
        std::vector<long long> last;
        bool ok = false;
        auto dfs = [&](auto&& self, int i) -> void {
            if (ok) return;
            if (i == n) {
                ok = true;
                return;
            }
            for (std::size_t c = 0; c < last.size(); ++c) {
                const bool fit = increasing ? less_by(last[c], v[i], strict) : less_by(v[i], last[c], strict);
                if (fit) {
                    const long long keep = last[c];
                    last[c] = v[i];
                    self(self, i + 1);
                    last[c] = keep;
                }
            }
            if ((int)last.size() < k) {
                last.push_back(v[i]);
                self(self, i + 1);
                last.pop_back();
            }
        };
        dfs(dfs, 0);
        if (ok) return k;
    }
    return n;
}

static std::vector<long long> rnd_values(std::mt19937_64& rng, int n) {
    std::vector<long long> v(n);
    const int t = (int)(rng() % 4);
    const long long range = 1 + (long long)(rng() % (t == 0 ? 3 : 2 * n + 2));
    for (int i = 0; i < n; ++i) {
        if (t == 3) v[i] = rnd_ll(rng);
        else v[i] = (long long)(rng() % range) - (t == 2 ? range / 2 : 0);
    }
    if (rng() % 8 == 0) std::sort(v.begin(), v.end());
    if (rng() % 8 == 0) std::sort(v.rbegin(), v.rend());
    return v;
}

static void check_restore(const std::vector<long long>& v, const std::vector<int>& r, int L, bool strict) {
    LIS_CHECK((int)r.size() == L);
    for (std::size_t k = 0; k < r.size(); ++k) {
        LIS_CHECK(r[k] >= 0 && r[k] < (int)v.size());
        if (k) LIS_CHECK(r[k - 1] < r[k] && less_by(v[r[k - 1]], v[r[k]], strict));
    }
}

static void test_small_exhaustive(std::mt19937_64& rng, int& cases) {
    // 全部の部分集合を調べる：LIS の長さ・数・入りうる位置・必ず入る位置
    for (int it = 0; it < 3000 * TEST_SCALE; ++it) {
        const int n = (int)(rng() % 13);
        const auto v = rnd_values(rng, n);
        const bool strict = rng() % 2;
        int L = 0;
        long long count = 0;
        std::vector<int> through(n, 0);
        for (int mask = 0; mask < (1 << n); ++mask) {
            int prev = -1, len = 0;
            bool ok = true;
            for (int i = 0; i < n && ok; ++i)
                if (mask >> i & 1) {
                    if (prev >= 0 && !less_by(v[prev], v[i], strict)) ok = false;
                    prev = i, ++len;
                }
            if (!ok) continue;
            if (len > L) L = len, count = 0, std::fill(through.begin(), through.end(), 0);
            if (len == L) {
                ++count;
                for (int i = 0; i < n; ++i)
                    if (mask >> i & 1) ++through[i];
            }
        }
        LIS_CHECK(lis_length(v, strict) == L);
        check_restore(v, lis_restore(v, strict), L, strict);
        LIS_CHECK(count_lis<long long>(v, strict) == count);
        LIS_CHECK(count_lis<Mod>(v, strict) == Mod(count));
        const LisMembers mb = lis_members(v, strict);
        LIS_CHECK(mb.length == L && (int)mb.in_some.size() == n && (int)mb.in_all.size() == n);
        for (int i = 0; i < n; ++i) LIS_CHECK(mb.in_some[i] == (through[i] > 0) && mb.in_all[i] == (L > 0 && through[i] == count));
        if (n <= 9) {
            LIS_CHECK(min_increasing_cover(v, strict) == bf_cover(v, strict, true));
            LIS_CHECK(min_decreasing_cover(v, strict) == bf_cover(v, strict, false));
        }
        ++cases;
    }
}

static void test_dp(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 1500 * TEST_SCALE; ++it) {
        const int n = (int)(rng() % (it % 20 == 0 ? 600 : 60));
        const auto v = rnd_values(rng, n);
        for (int st = 0; st < 2; ++st) {
            const bool strict = st == 1;
            const Bf b = bf_dp(v, strict);
            LIS_CHECK(lis_length(v, strict) == b.L);
            LIS_CHECK(lis_ending_at(v, strict) == b.e && lis_starting_at(v, strict) == b.s);
            check_restore(v, lis_restore(v, strict), b.L, strict);
            LIS_CHECK(count_lis<Mod>(v, strict) == b.total);
            // 減少：大小を逆にした LIS
            std::vector<long long> w(v);
            for (auto& x : w) x = ~x;
            LIS_CHECK(lds_length(v, strict) == bf_dp(w, strict).L);
            LIS_CHECK(min_increasing_cover(v, strict) == lds_length(v, !strict) && min_decreasing_cover(v, strict) == lis_length(v, !strict));
            const LisMembers mb = lis_members(v, strict);
            for (int i = 0; i < n; ++i) {
                const bool some = b.e[i] + b.s[i] - 1 == b.L;
                // すべての LIS に入る ⇔ i を通る LIS の数（ce[i] × cs[i]）が全体の数と同じ（mod で比べる）
                const Mod through((long long)((unsigned long long)b.ce[i].v * b.cs[i].v % 998244353ULL));
                LIS_CHECK(mb.in_some[i] == some);
                LIS_CHECK(mb.in_all[i] == (some && through == b.total));
            }
        }
        // 文字列でも同じ
        const std::string s = rnd_str_shaped(rng, (int)(rng() % 40), 1 + (int)(rng() % 26));
        const std::vector<long long> sv(s.begin(), s.end());
        LIS_CHECK(lis_length(s) == lis_length(sv) && lis_length(s, false) == lis_length(sv, false) && lis_restore(s) == lis_restore(sv));
        ++cases;
    }
}

static void test_2d(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 1500 * TEST_SCALE; ++it) {
        const int n = (int)(rng() % 60);
        const long long range = 1 + (long long)(rng() % 12);
        std::vector<std::pair<long long, int>> ps(n);
        for (auto& p : ps) p = {(long long)(rng() % range) - 3, (int)(rng() % range)};
        // 愚直：x で並べて O(N^2)
        std::vector<int> ord(n);
        for (int i = 0; i < n; ++i) ord[i] = i;
        std::sort(ord.begin(), ord.end(), [&](int i, int j) { return ps[i].first < ps[j].first; });
        std::vector<int> dp(n, 1);
        int L = 0;
        for (int a = 0; a < n; ++a) {
            for (int b = 0; b < a; ++b)
                if (ps[ord[b]].first < ps[ord[a]].first && ps[ord[b]].second < ps[ord[a]].second) dp[a] = std::max(dp[a], dp[b] + 1);
            L = std::max(L, dp[a]);
        }
        const auto r = lis_2d(ps);
        LIS_CHECK((int)r.size() == L);
        for (std::size_t k = 1; k < r.size(); ++k) LIS_CHECK(ps[r[k - 1]].first < ps[r[k]].first && ps[r[k - 1]].second < ps[r[k]].second);
        ++cases;
    }
}

static void random_test() {
    std::mt19937_64 rng(20261012);
    int cases = 0;
    test_small_exhaustive(rng, cases);
    test_dp(rng, cases);
    test_2d(rng, cases);
    // 途中の数が long long を超えても、答えが収まれば止まらない
    {
        std::vector<long long> v;
        for (int k = 0; k <= 62; ++k) v.push_back(1000 + 2 * k + 1), v.push_back(1000 + 2 * k);
        for (int k = 0; k < 64; ++k) v.push_back(k);
        LIS_CHECK(lis_length(v) == 64 && count_lis<long long>(v) == 1);
        ++cases;
    }
    // lis_2d は < だけで比べる（double・unsigned long long）
    {
        LIS_CHECK(lis_2d(std::vector<std::pair<double, double>>{{0, 0.1}, {1, 0.5}}).size() == 2);
        LIS_CHECK(lis_2d(std::vector<std::pair<long long, unsigned long long>>{{0, 1}, {1, (1ULL << 63) + 5}}).size() == 2);
        ++cases;
    }
    std::cout << "[lis] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[lis] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    auto ms = [](std::chrono::steady_clock::time_point a) {
        return (long long)std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count();
    };
    const int N = 1000000;
    std::vector<int> a(N);
    for (auto& x : a) x = (int)(rng() % 1000000000);
    auto t0 = std::chrono::steady_clock::now();
    long long x = lis_length(a);
    const long long t_len = ms(t0);
    t0 = std::chrono::steady_clock::now();
    x += (long long)lis_restore(a).size();
    const long long t_res = ms(t0);
    t0 = std::chrono::steady_clock::now();
    const LisMembers mb = lis_members(a);
    x += mb.in_some[N / 2];
    const long long t_mem = ms(t0);
    t0 = std::chrono::steady_clock::now();
    x += count_lis<Mod>(a).v;
    const long long t_cnt = ms(t0);
    std::vector<std::pair<int, int>> ps(N);
    for (auto& p : ps) p = {(int)(rng() % 1000000), (int)(rng() % 1000000)};
    t0 = std::chrono::steady_clock::now();
    x += (long long)lis_2d(ps).size();
    const long long t_2d = ms(t0);
    std::cout << "[lis] benchmark: N = 10^6  lis_length " << t_len << " ms, lis_restore " << t_res << " ms, lis_members " << t_mem << " ms, count_lis "
              << t_cnt << " ms, lis_2d " << t_2d << " ms  (check " << x % 7 << ")\n";
}

}  // namespace lis_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace lis_demo {

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
#define SHOW(expr, comment) lis_demo::show(#expr, (expr), comment)
#define RUN(...)                   \
    lis_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("最長増加部分列（lis） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  lis_length / lis_restore / lis_ending_at / lis_starting_at   長さ・1 つ（添字）・各位置で終わる / 始まる長さ\n"
              << "  lis_members / count_lis<T>                                 入りうる・必ず入る位置 / LIS の数\n"
              << "  lds_length / min_increasing_cover / min_decreasing_cover    減少・増加列に分ける最小の数（Dilworth）\n"
              << "  lis_2d                                                     (x, y) が両方とも増える最長の列（入れ子の箱）\n"
              << "  strict = true（既定）で狭義（<）、false で広義（≤）。列は vector<整数> でも string でもよい\n";

    const std::vector<int> a = {3, 1, 4, 1, 5, 9, 2, 6, 5, 3, 5};
    code("const vector<int> a = {3, 1, 4, 1, 5, 9, 2, 6, 5, 3, 5};");

    api("【1】 lis_length(a, strict) / lis_restore(a, strict)", "int / vector<int>（添字。昇順）", "最長増加部分列の長さと、その 1 つ", "O(N log N)");
    SHOW(lis_length(a), "3, 4, 5, 6 など（狭義）");
    SHOW(lis_length(a, false), "1, 1, 2, 3, 5 など（広義 = 減らない）");
    SHOW(lis_restore(a), "a[3], a[6], a[9], a[10] = 1, 2, 3, 5（いくつかあるうちの 1 つ）");
    SHOW(lis_length(std::string("abracadabra")), "文字列でも（a < b < c < d < r）");

    api("【2】 lis_ending_at(a, strict) / lis_starting_at(a, strict)", "vector<int>", "各位置で終わる / 始まる最長の増加部分列の長さ（両方を足して 1 を引くと、その位置を通る最長）", "O(N log N)");
    SHOW(lis_ending_at(a), "");
    SHOW(lis_starting_at(a), "");

    api("【3】 lis_members(a, strict)", "LisMembers{length, in_some, in_all}", "各位置が「どれかの LIS に入る」「すべての LIS に入る」か", "O(N log N)");
    {
        const LisMembers mb = lis_members(a);
        SHOW(mb.in_some, "ABC354 F（Useless for LIS）はこれが true の位置");
        SHOW(mb.in_all, "");
        const std::vector<int> c = {2, 5, 3, 4, 3, 4};
        const LisMembers m2 = lis_members(c);
        std::vector<int> idx;
        for (int i = 0; i < (int)c.size(); ++i)
            if (m2.in_some[i]) idx.push_back(i + 1);
        note("ABC354 F の入力例 2 の 1 つ目（2 5 3 4 3 4）→ " + to_s((int)idx.size()) + " 個: " + join(idx));
    }

    api("【4】 count_lis<T>(a, strict)", "T", "LIS の数（添字の列として。同じ値の並びでも位置が違えば別）。T は ModInt998・long long（溢れたら止める）など", "O(N log N)");
    SHOW(count_lis<long long>(a), "");
    SHOW(count_lis<long long>(std::vector<int>{2, 1, 4, 3, 6, 5}), "2 × 2 × 2");

    api("【5】 lds_length(a, strict) / min_increasing_cover(a, strict) / min_decreasing_cover(a, strict)", "int",
        "最長減少部分列 / 増加（減少）部分列に分ける最小の数（Dilworth：狭義増加に分ける数 = 最長の「減らない」部分列）", "O(N log N)");
    SHOW(lds_length(a), "9, 6, 5, 3 など");
    SHOW(min_increasing_cover(std::vector<int>{2, 1, 4, 5, 3}), "ABC134 E（Sequence Decomposing）の入力例 1");
    SHOW(min_increasing_cover(std::vector<int>{0, 0, 0, 0}), "入力例 2：同じ値は同じ色にできない");

    api("【6】 lis_2d(ps)", "vector<int>（添字。x の小さい順）", "(x, y) が両方とも狭義に増える最長の列（x の昇順・同じ x は y の降順に並べて y の LIS）", "O(N log N)");
    {
        const std::vector<std::pair<int, int>> box = {{8, 8}, {5, 3}, {2, 2}, {4, 2}, {2, 1}};
        code("const vector<pair<int, int>> box = {{8, 8}, {5, 3}, {2, 2}, {4, 2}, {2, 1}};   // ABC038 D（プレゼント）の入力例 4");
        SHOW(lis_2d(box).size(), "入れ子にできる箱の最大の数");
        SHOW(lis_2d(box), "(2, 1) → (4, 2) → (5, 3) → (8, 8)");
    }

    head("問題での使い方：ABC360 G（Suitable Edit for LIS）");
    note("1 か所を好きな整数に変えて LIS を最大に → e = lis_ending_at(A)、s = lis_starting_at(A) を作り、");
    note("j < k、k - j ≥ 2、A[k] - A[j] ≥ 2 の組で e[j] + 1 + s[k] の最大（間の 1 つを A[j] + 1 に変える）。両端の外も番兵で。");
    note("本番は j を値の BIT（最大値）に入れて O(N log N)。ここでは入力例を O(N^2) で:");
    {
        for (const auto& A : {std::vector<long long>{3, 2, 2, 4}, std::vector<long long>{4, 5, 3, 6, 7}}) {
            const int n = (int)A.size();
            const auto e = lis_ending_at(A), s = lis_starting_at(A);
            int best = 0;
            for (int x = 0; x < n; ++x) {  // x を変える。左は j < x、右は k > x
                int left = 0, right = 0;
                for (int j = 0; j < x; ++j) left = std::max(left, e[j]);
                for (int k = x + 1; k < n; ++k) right = std::max(right, s[k]);
                best = std::max({best, left + 1, 1 + right});
                for (int j = 0; j < x; ++j)
                    for (int k = x + 1; k < n; ++k)
                        if (A[k] - A[j] >= 2) best = std::max(best, e[j] + 1 + s[k]);
            }
            note("A = " + to_s(A) + " → " + to_s(best));
        }
    }
}

}  // namespace lis_demo

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
    lis_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    lis_test::benchmark();    // 速度計測
    lis_demo::demo();         // 使い方デモ
    return 0;
}
