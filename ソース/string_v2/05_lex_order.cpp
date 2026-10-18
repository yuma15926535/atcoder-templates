// =============================================================================
//  辞書順（つなげて最小・独自の文字の順・最小 / 最大の回転・列名 a, …, z, aa・辞書の番号・同じ長さの次の文字列）   --- 文字列第二世代
// -----------------------------------------------------------------------------
//  分類       : 05 辞書順・回転
//  依存       : 00_str_base（同梱済み）
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）。__int128 を使う（GCC / Clang）
//  文字コード : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  辞書順は std::string と同じ（先頭から比べ、片方が接頭辞なら短い方が小さい。char は 0〜255 として比べる）
//  文字の種類（列名・辞書の番号・次の文字列）は Alphabet で指定する（既定は英小文字）
//
//  機能一覧                                                                       計算量
//    concat_order(vs) / min_concat(vs)     全部をつなげて辞書順最小になる並べ方（添字）/ その列（a + b < b + a）   O(N log N × 最長の長さ)
//    lex_less_by(order) / lex_less_by(A)   独自の文字の順の比較関数 LexLess（sort に渡す。.compare(a, b) は -1 / 0 / 1）
//    min_rotation(s) / max_rotation(s)     辞書順最小 / 最大の回転 s[k..] + s[..k] の k（同じなら最小の k）           O(N)
//    bijective_to_string(n, A)             列名：1 → a、26 → z、27 → aa（0 なら空）                                 O(log n)
//    bijective_from_string(s, A)           その逆（a → 1）                                                         O(|s|)
//    dict_rank(s, max_len, A)              長さ 1〜max_len の全部の文字列を辞書順に並べたときの s の番号（0 から）   O(|s| + max_len)
//    dict_unrank(k, max_len, A)            その逆（番号 k の文字列）                                                 O(max_len)
//    next_string(s, A)                     同じ長さで辞書順の次の文字列（最後なら空）                                O(N)
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    lex_test::random_test()  … 全部の並べ方・全部の文字列を並べる愚直解と全機能を突き合わせるランダムテスト
//    lex_test::benchmark()    … 速度計測
//    lex_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
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
#ifndef GG_LEX_ORDER
#define GG_LEX_ORDER
// ---- 辞書順（文字列第二世代）  つなげて最小・独自の文字の順・最小の回転・列名（a, …, z, aa）・辞書の番号・同じ長さでの番号・次の文字列 ----
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
template <class T, class Seq>
T lex_rank_same_length(const Seq& s0, const Alphabet& A = Alphabet::lower()) {  // A の文字からなる長さ |s| の文字列で s より小さいものの数（= s を σ 進数で読んだ値）  O(N)
    const auto& s = str_internal::seq(s0);
    T r = T(0);
    for (std::size_t i = 0; i < s.size(); ++i)
        r = str_internal::checked_add(str_internal::checked_mul(r, T(A.size()), "lex_rank_same_length"), T(A.id(s[i])), "lex_rank_same_length");
    return r;
}
#endif  // GG_LEX_ORDER
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト
// =============================================================================
namespace lex_test {

#define LEX_CHECK(...)                                                                          \
    do {                                                                                        \
        if (!(__VA_ARGS__)) {                                                                   \
            std::cerr << "[lex_order] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
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

static void test_concat(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 3000 * TEST_SCALE; ++it) {
        const int k = (int)(rng() % 6);
        std::vector<std::string> vs(k);
        for (auto& s : vs) s = rnd_str_shaped(rng, (int)(rng() % 4), 1 + (int)(rng() % 2));
        std::vector<int> p(k);
        std::iota(p.begin(), p.end(), 0);
        std::string best;
        bool first = true;
        do {
            std::string t;
            for (int i : p) t += vs[i];
            if (first || t < best) best = t, first = false;
        } while (std::next_permutation(p.begin(), p.end()));
        LEX_CHECK(min_concat(vs) == best);
        const auto ord = concat_order(vs);
        std::string t;
        for (int i : ord) t += vs[i];
        LEX_CHECK(t == best);
        ++cases;
    }
}

static void test_order(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 2000 * TEST_SCALE; ++it) {
        std::string ord = "abcdefghijklmnopqrstuvwxyz";
        std::shuffle(ord.begin(), ord.end(), rng);
        const LexLess lt = lex_less_by(ord);
        std::vector<int> pos(256);
        for (int i = 0; i < 26; ++i) pos[(unsigned char)ord[i]] = i;
        for (int t = 0; t < 20; ++t) {
            const std::string a = rnd_str(rng, (int)(rng() % 4), 3), b = rnd_str(rng, (int)(rng() % 4), 3);
            std::vector<int> x, y;
            for (char c : a) x.push_back(pos[(unsigned char)c]);
            for (char c : b) y.push_back(pos[(unsigned char)c]);
            LEX_CHECK(lt(a, b) == (x < y) && lt.compare(a, b) == (x < y ? -1 : (x == y ? 0 : 1)));
            ++cases;
        }
    }
}

static void test_rotation(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 5000 * TEST_SCALE; ++it) {
        const std::string s = rnd_str_shaped(rng, (int)(rng() % 14), 1 + (int)(rng() % 3));
        const int n = (int)s.size();
        int mn = 0, mx = 0;
        std::string bmn = s, bmx = s;
        for (int k = 0; k < n; ++k) {
            const std::string r = s.substr(k) + s.substr(0, k);
            if (r < bmn) bmn = r, mn = k;
            if (r > bmx) bmx = r, mx = k;
        }
        LEX_CHECK(min_rotation(s) == mn && max_rotation(s) == mx);
        LEX_CHECK(min_rotation(std::vector<int>(s.begin(), s.end())) == mn);
        ++cases;
    }
}

static void test_names(std::mt19937_64& rng, int& cases) {
    // 列名：1 から順に作ったものと
    std::string cur;
    for (long long n = 1; n <= 20000; ++n) {
        // cur を「次」に（z の次は aa）
        int i = (int)cur.size() - 1;
        while (i >= 0 && cur[i] == 'z') cur[i--] = 'a';
        if (i < 0) cur.insert(cur.begin(), 'a');
        else ++cur[i];
        LEX_CHECK(bijective_to_string(n) == cur && bijective_from_string(cur) == n);
        ++cases;
    }
    LEX_CHECK(bijective_to_string(0).empty() && bijective_to_string(28, Alphabet::upper()) == "AB");
    for (int it = 0; it < 1000 * TEST_SCALE; ++it) {
        const long long n = rnd_ll(rng) & std::numeric_limits<long long>::max();
        LEX_CHECK(bijective_from_string(bijective_to_string(n)) == n);
        const Alphabet A = Alphabet::of("xy");
        LEX_CHECK(bijective_from_string(bijective_to_string(n, A), A) == n);
    }
    // 辞書の番号：全部の文字列を並べたものと
    for (int sg = 1; sg <= 3; ++sg)
        for (int L = 1; L <= 6; ++L) {
            std::string chars = std::string("abc").substr(0, sg);
            const Alphabet A = Alphabet::of(chars);
            std::vector<std::string> all;
            std::vector<std::string> layer = {""};
            for (int len = 1; len <= L; ++len) {
                std::vector<std::string> nl;
                for (const auto& w : layer)
                    for (char c : chars) nl.push_back(w + c);
                all.insert(all.end(), nl.begin(), nl.end());
                layer = nl;
            }
            std::sort(all.begin(), all.end());
            for (long long k = 0; k < (long long)all.size(); ++k) {
                LEX_CHECK(dict_rank(all[k], L, A) == k && dict_unrank(k, L, A) == all[k]);
                ++cases;
            }
            // 同じ長さの次
            for (const auto& w : layer) {
                auto it2 = std::upper_bound(layer.begin(), layer.end(), w);
                const std::string nx = it2 == layer.end() ? "" : *it2;
                LEX_CHECK(next_string(w, A) == nx);
            }
        }
    // 大きな max_len
    LEX_CHECK(dict_rank("aaa", 1000000000000LL, Alphabet::of("a")) == 2 && dict_rank("aab", 3) == 3 && dict_unrank(0, 1000000000000000000LL) == "a");
    LEX_CHECK(dict_rank("ab", 40, Alphabet::of("ab")) == (1LL << 39));  // "a" と、"aa" で始まる 2^39 - 1 個
    LEX_CHECK(dict_unrank(1LL << 39, 40, Alphabet::of("ab")) == "ab");
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
static void test_final_additions(std::mt19937_64& rng, int& cases) {  // lex_rank_same_length
    for (int it = 0; it < 1500 * TEST_SCALE; ++it) {
        const int sg = 1 + (int)(rng() % 3), n = (int)(rng() % 7);
        std::string letters = std::string("xyz").substr(0, sg);
        std::shuffle(letters.begin(), letters.end(), rng);
        const Alphabet A = Alphabet::of(letters);
        std::string s(n, ' ');
        for (char& c : s) c = letters[rng() % sg];
        const std::vector<int> sid = to_ids(s, A);
        long long cnt = 0;
        for (const auto& v : all_id_strings(n, sg)) cnt += v < sid;
        LEX_CHECK(lex_rank_same_length<long long>(s, A) == cnt && lex_rank_same_length<double>(s, A) == (double)cnt);
        ++cases;
    }
    LEX_CHECK(lex_rank_same_length<long long>(std::vector<int>{1, 0, 1}, Alphabet::ints(2)) == 5);
}
static void random_test() {
    std::mt19937_64 rng(20261012);
    int cases = 0;
    test_concat(rng, cases);
    test_order(rng, cases);
    test_rotation(rng, cases);
    test_names(rng, cases);
    test_final_additions(rng, cases);
    std::cout << "[lex_order] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[lex_order] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    auto ms = [](std::chrono::steady_clock::time_point a) {
        return (long long)std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count();
    };
    const int n = 10000000;
    const std::string s = rnd_str_shaped(rng, n, 2), a(n, 'a');
    auto t0 = std::chrono::steady_clock::now();
    const long long r = min_rotation(s) + max_rotation(s) + min_rotation(a);
    const long long t_rot = ms(t0);
    std::vector<std::string> vs(200000);
    for (auto& x : vs) x = rnd_str(rng, 1 + (int)(rng() % 10), 3);
    t0 = std::chrono::steady_clock::now();
    const auto ord = concat_order(vs);
    const long long t_cat = ms(t0);
    std::vector<std::string> ws(200000);
    for (auto& x : ws) x = rnd_str(rng, 10, 26);
    t0 = std::chrono::steady_clock::now();
    std::sort(ws.begin(), ws.end(), lex_less_by("qwertyuiopasdfghjklzxcvbnm"));
    const long long t_sort = ms(t0);
    std::cout << "[lex_order] benchmark: min_rotation + max_rotation + min_rotation（N = 10^7）" << t_rot << " ms, concat_order（2×10^5 個）" << t_cat
              << " ms, lex_less_by で sort（2×10^5 個）" << t_sort << " ms  (check " << (r + ord[5] + ws[3][0]) % 7 << ")\n";
}

}  // namespace lex_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace lex_demo {

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
#define SHOW(expr, comment) lex_demo::show(#expr, (expr), comment)
#define RUN(...)                   \
    lex_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("辞書順（lex_order） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  concat_order / min_concat                 つなげて辞書順最小（a + b < b + a の順）\n"
              << "  lex_less_by(order) / LexLess::compare    独自の文字の順の比較（sort に渡す）\n"
              << "  min_rotation / max_rotation              辞書順最小・最大の回転の開始位置\n"
              << "  bijective_to_string / bijective_from_string   列名（a, …, z, aa, …）と番号\n"
              << "  dict_rank / dict_unrank                  長さ max_len 以下の全部の文字列での辞書順の番号\n"
              << "  next_string                              同じ長さで辞書順の次\n"
              << "  lex_rank_same_length<T>(s, A)            同じ長さの文字列の中での辞書順の番号（mod でも）\n";

    api("【1】 concat_order(vs) / min_concat(vs)", "vector<int> / 列",
        "全部をつなげて辞書順最小になる並べ方（a + b < b + a で並べる。単に辞書順に並べるのとは違う）", "O(N log N × 長さ)");
    RUN(const std::vector<std::string> vs = {"b", "ba", "bab"});
    SHOW(min_concat(vs), "辞書順の b, ba, bab → \"bbabab\" より小さい");
    SHOW(concat_order(vs), "");
    note("ABC225 F（String Cards）: concat_order で並べてから、後ろから「K 枚選ぶ」DP（dp[i][j] = min(dp[i+1][j], S_i + dp[i+1][j-1])）。");

    api("【2】 lex_less_by(order) / LexLess::compare(a, b)", "比較関数",
        "order に並べた文字の順での辞書順（order にない文字は止まる）。sort にそのまま渡せる", "O(長さ)");
    {
        std::vector<std::string> names = {"abx", "bzz", "bzy", "caa"};
        std::sort(names.begin(), names.end(), lex_less_by("bacdefghijklmnopqrstuvwxzy"));
        note("ABC219 C（Neo-lexicographic Ordering）の入力例 1：sort(names, lex_less_by(X))");
        std::cout << "    → " << join(names) << "\n";
    }
    SHOW(lex_less_by("cba").compare("ab", "ba"), "c < b < a なので \"ab\" > \"ba\"");

    api("【3】 min_rotation(s) / max_rotation(s)", "int",
        "辞書順最小 / 最大の回転 s[k..] + s[..k] の k（同じ回転が何度も出るなら最小の k）", "O(N)");
    RUN(const std::string S = "cabca");
    SHOW(min_rotation(S), "\"acabc\"");
    SHOW(max_rotation(S), "\"cabca\"");
    {
        // ABC223 B: 回転の中で辞書順最小と最大
        const std::string T = "aaba";
        const int a = min_rotation(T), b = max_rotation(T);
        note("ABC223 B（String Shifting）: min_rotation と max_rotation の位置で回したもの");
        std::cout << "    T = aaba → " << T.substr(a) + T.substr(0, a) << " / " << T.substr(b) + T.substr(0, b) << "\n";
    }

    api("【4】 bijective_to_string(n, A) / bijective_from_string(s, A)", "string / long long",
        "列名（全単射の記数法）：1 → a、26 → z、27 → aa、28 → ab。0 は空。A で文字を変えられる", "O(log n)");
    SHOW(bijective_to_string(2), "ABC171 C（One Quadrillion and One Dalmatians）の入力例 1");
    SHOW(bijective_to_string(27), "");
    SHOW(bijective_to_string(123456789), "ABC171 C の入力例 3");
    SHOW(bijective_from_string("aa"), "ABC285 C");
    SHOW(bijective_from_string("BRUTMHYHIIZP", Alphabet::upper()), "ABC285 C の入力例 3");

    api("【5】 dict_rank(s, max_len, A) / dict_unrank(k, max_len, A)", "long long / string",
        "長さ 1〜max_len の全部の文字列を辞書順に並べたときの番号（0 から）とその逆。a, aa, aaa, …, ab, … の順", "O(max_len)");
    SHOW(dict_rank("a", 3), "最初");
    SHOW(dict_rank("ab", 3), "a, aa, aaa, …, aaz, ab");
    SHOW(dict_unrank(28, 3), "");
    SHOW(dict_rank("ba", 2, Alphabet::of("ab")), "a, aa, ab, b, ba");

    api("【6】 next_string(s, A)", "string", "同じ長さで辞書順の次（\"az\" → \"ba\"。最後の \"zz…\" なら空）", "O(N)");
    SHOW(next_string("az"), "");
    SHOW(next_string("zz"), "ない → 空");
    SHOW(next_string("011", Alphabet::of("01")), "2 進数の +1");

    api("【7】 lex_rank_same_length<T>(s, A)", "T（mint・long long。long long で溢れたら止める）",
        "A の文字からなる長さ |s| の文字列で s より小さいものの数（s を σ 進数で読んだ値。N = 10^6 でも mint で）", "O(N)");
    SHOW(lex_rank_same_length<long long>("ba"), "aa〜az の 26 個");
    SHOW(lex_rank_same_length<long long>(std::vector<int>{1, 0, 1}, Alphabet::ints(2)), "2 進数の 101");
    note("長さ max_len 以下の全部の文字列での順位は dict_rank。同じ文字の個数が決まった並べ替えの中での順位・k 番目は math_v2 04_permutation の multiset_perm_rank / multiset_perm_unrank（ABC202 D・ABC215 C）。");
}

}  // namespace lex_demo

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
    lex_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    lex_test::benchmark();    // 速度計測
    lex_demo::demo();         // 使い方デモ
    return 0;
}
