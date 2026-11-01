// =============================================================================
//  1 文字の書き換えがある文字列（区間のハッシュ・比較・LCP・回文判定・文字の数・並んでいるか・k 番目の文字）   --- 文字列第二世代
// -----------------------------------------------------------------------------
//  分類       : 09 操作・更新のある文字列
//  依存       : 00_str_base、01_rolling_hash（同梱済み）
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）
//  文字コード : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  ハッシュは 01_rolling_hash と同じ基数・同じ値（RollingHash・hash_of で作ったほかの文字列とも比べられる）。中はフェニック木
//  文字の数・並び（count・counts・is_sorted・kth_char）の表は初めて使うときに作る（メモリ (N+1)σ 個の int。そのときに A の範囲を確かめる）
//  列は string でも vector<int>・vector<long long> でもよい。範囲は [l, r)。「L 文字目から R 文字目」は (L - 1, R)
//
//  機能一覧                                                                       計算量
//    MutableString ms(s, A)                作る（A は文字の種類。既定は英小文字）                          O(N)
//      ms.set(i, c) / get(i)               i 文字目を c に / i 文字目の値（char は 0〜255）                 O(log N) / O(1)
//      ms.hash(l, r) / reversed_hash(l, r) s[l, r) / 反転した s[l, r) のハッシュ（StrHash）                  O(log N)
//      ms.same(l1, r1, l2, r2)             s[l1, r1) == s[l2, r2) か                                    O(log N)
//      ms.is_palindrome(l, r)              s[l, r) が回文か（ABC331 F）                                   O(log N)
//      ms.lcp(i, j) / compare(l1, r1, l2, r2)   最長共通接頭辞 / 辞書順（-1 / 0 / 1）                      O(log^2 N)
//      ms.count(l, r, c) / counts(l, r)    s[l, r) にある c の数 / 番号ごとの数                            O(log N) / O(σ log N)
//      ms.is_sorted(l, r)                  s[l, r) が広義単調増加（A の番号の順）か（ABC285 F）              O(log N)
//      ms.kth_char(l, r, k)                s[l, r) を並べ替えた k 番目の値                                 O(σ log N)
//      ms.size() / values() / to_string()
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    ms_test::random_test()  … 書き換えのたびに作り直した文字列と全機能を突き合わせるランダムテスト
//    ms_test::benchmark()    … 速度計測
//    ms_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
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
#ifndef GG_ROLLING_HASH
#define GG_ROLLING_HASH
// ---- ローリングハッシュ（文字列第二世代）  mod 2^61 - 1。基数は実行ごとの乱数で、すべての文字列で共通 ----
// 別の文字列の部分文字列どうし・つなげたものも比べられる。累乗の表は基数を変えずに伸ばす
// 1 回の比較で衝突する確率は約 N / 2^61（N は長さ）。k 文字まで書き換えてよい回文の長さ（center_palindrome_lengths_with_changes）もここ
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

// ---- k 文字まで書き換えてよい回文：中心 c（0〜2N-2。偶数 c は s[c/2] が中心の奇数長、奇数 c は間が中心の偶数長）ごとの最長の長さ ----
//   中心から広げ、食い違いを k 回まで飛ばす（ハッシュの LCP を k + 1 回）。k = 0 なら 04_palindrome の center_palindrome_lengths と同じ（ABC468 D）
//   中心 c で書き換え k 回以内の回文は、長さ L[c] 以下で同じ偶奇のもの全部  O(N k log N)（k が小さいとき用。N = 2×10^5 なら k は 10 程度まで）
template <class Seq>
std::vector<int> center_palindrome_lengths_with_changes(const Seq& s0, int k) {
    if (k < 0) str_internal::fail("center_palindrome_lengths_with_changes", "k は 0 以上");
    const auto& s = str_internal::seq(s0);
    const int n = (int)s.size();
    if (n == 0) return {};
    const RollingHash F(s), B(reversed(s));
    std::vector<int> res(2 * n - 1);
    for (int c = 0; c < 2 * n - 1; ++c) {
        int l = c % 2 == 0 ? c / 2 - 1 : (c - 1) / 2, r = c % 2 == 0 ? c / 2 + 1 : (c + 1) / 2;  // 次に比べる左右の位置
        int len = c % 2 == 0 ? 1 : 0;
        for (int used = 0;; ++used) {
            if (l < 0 || r >= n) break;
            const int e = F.lcp(r, B, n - 1 - l);  // s[r..) と s[..l] を逆から読んだものの一致
            l -= e, r += e, len += 2 * e;
            if (l < 0 || r >= n || used == k) break;
            --l, ++r, len += 2;  // 食い違いを 1 つ書き換える
        }
        res[c] = len;
    }
    return res;
}

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
#ifndef GG_MUTABLE_STRING
#define GG_MUTABLE_STRING
// ---- 1 文字の書き換えがある文字列（文字列第二世代） ----
// ハッシュは 01_rolling_hash と同じ基数（RollingHash・hash_of で作ったほかの文字列とも比べられる）。フェニック木で O(log N)
// 文字の数・並び（count・counts・is_sorted・kth_char）の表は初めて使うときに作る（メモリ (N+1)σ 個の int。そのときに A の範囲を確かめる）
struct MutableString {
    MutableString() = default;
    template <class Seq>
    explicit MutableString(const Seq& s0, const Alphabet& A = Alphabet::lower()) : A_(A) {
        using C = str_internal::RHCore;
        const auto& s = str_internal::seq(s0);
        n_ = (int)s.size();
        C::ensure((std::size_t)n_ + 1);
        const unsigned long long inv = C::pow_mod(C::base1, C::MOD - 2);
        ipw_.assign(n_ + 1, 1);
        for (int i = 1; i <= n_; ++i) ipw_[i] = C::mul(ipw_[i - 1], inv);
        val_.resize(n_);
        fw_.assign(n_ + 1, 0), bw_.assign(n_ + 1, 0);
        for (int i = 0; i < n_; ++i) {
            val_[i] = str_internal::code_of(s[i]);
            const unsigned long long h = C::value(val_[i]);
            fw_[i + 1] = C::mul(h, ipw_[i]), bw_[i + 1] = C::mul(h, C::pw1[i]);
        }
        for (int i = 1; i <= n_; ++i) {  // フェニック木を O(N) で作る
            const int j = i + (i & -i);
            if (j <= n_) fw_[j] = C::add(fw_[j], fw_[i]), bw_[j] = C::add(bw_[j], bw_[i]);
        }
    }
    int size() const { return n_; }
    long long get(int i) const {  // i 文字目の値（char は 0〜255）  O(1)
        check_pos(i, "MutableString::get");
        return val_[i];
    }
    template <class T>
    void set(int i, T c) {  // i 文字目を c に  O(log N)
        using C = str_internal::RHCore;
        check_pos(i, "MutableString::set");
        const long long nv = str_internal::code_of(c), ov = val_[i];
        if (nv == ov) return;
        if (counts_built_) {
            const int oc = A_.id(ov), nc = A_.id(nv);
            bit_add(cnt_, oc, i, -1), bit_add(cnt_, nc, i, +1);
        }
        val_[i] = nv;
        if (desc_built_) {  // 前後との「下がる所」を付け直す
            for (int k = std::max(0, i - 1); k <= std::min(n_ - 2, i); ++k) {
                const int d = A_.id(val_[k]) > A_.id(val_[k + 1]) ? 1 : 0;
                if (d != desc_flag_[k]) bit_add_desc(k, d - desc_flag_[k]), desc_flag_[k] = (char)d;
            }
        }
        const unsigned long long diff = C::sub(C::value(nv), C::value(ov));
        for (int k = i + 1; k <= n_; k += k & -k) fw_[k] = C::add(fw_[k], C::mul(diff, ipw_[i])), bw_[k] = C::add(bw_[k], C::mul(diff, C::pw1[i]));
    }
    StrHash hash(int l, int r) const {  // s[l, r) のハッシュ（RollingHash::get と同じ値）  O(log N)
        str_internal::check_range("MutableString::hash", l, r, n_);
        return {raw(l, r), (long long)(r - l)};
    }
    StrHash reversed_hash(int l, int r) const {  // s[l, r) を反転した文字列のハッシュ  O(log N)
        str_internal::check_range("MutableString::reversed_hash", l, r, n_);
        return {raw_rev(l, r), (long long)(r - l)};
    }
    bool same(int l1, int r1, int l2, int r2) const {  // s[l1, r1) == s[l2, r2) か  O(log N)
        str_internal::check_range("MutableString::same", l1, r1, n_), str_internal::check_range("MutableString::same", l2, r2, n_);
        return r1 - l1 == r2 - l2 && hash(l1, r1) == hash(l2, r2);
    }
    bool is_palindrome(int l, int r) const {  // s[l, r) が回文か（ABC331 F）  O(log N)
        str_internal::check_range("MutableString::is_palindrome", l, r, n_);
        return raw(l, r) == raw_rev(l, r);
    }
    int lcp(int i, int j) const {  // s[i..) と s[j..) の最長共通接頭辞の長さ  O(log^2 N)
        if (i < 0 || i > n_ || j < 0 || j > n_) str_internal::fail("MutableString::lcp", "位置は 0〜N");
        int lo = 0, hi = std::min(n_ - i, n_ - j);
        while (lo < hi) {
            const int mid = (lo + hi + 1) >> 1;
            if (raw(i, i + mid) == raw(j, j + mid)) lo = mid;
            else hi = mid - 1;
        }
        return lo;
    }
    int compare(int l1, int r1, int l2, int r2) const {  // s[l1, r1) と s[l2, r2) の辞書順（値の大小で。-1 / 0 / 1）  O(log^2 N)
        str_internal::check_range("MutableString::compare", l1, r1, n_), str_internal::check_range("MutableString::compare", l2, r2, n_);
        const int m = std::min(r1 - l1, r2 - l2);
        int lo = 0, hi = m;
        while (lo < hi) {
            const int mid = (lo + hi + 1) >> 1;
            if (raw(l1, l1 + mid) == raw(l2, l2 + mid)) lo = mid;
            else hi = mid - 1;
        }
        if (lo == m) return (r1 - l1) == (r2 - l2) ? 0 : ((r1 - l1) < (r2 - l2) ? -1 : 1);
        return val_[l1 + lo] < val_[l2 + lo] ? -1 : 1;
    }

    // ---- 文字の数・並び（Alphabet の番号の順） ----
    template <class T>
    long long count(int l, int r, T c) const {  // s[l, r) にある c の数  O(log N)
        str_internal::check_range("MutableString::count", l, r, n_);
        ensure_counts();
        const int id = A_.id(c);
        return bit_sum(cnt_, id, r) - bit_sum(cnt_, id, l);
    }
    std::vector<long long> counts(int l, int r) const {  // s[l, r) の番号ごとの数（長さ σ）  O(σ log N)
        str_internal::check_range("MutableString::counts", l, r, n_);
        ensure_counts();
        std::vector<long long> res(A_.size());
        for (int c = 0; c < A_.size(); ++c) res[c] = bit_sum(cnt_, c, r) - bit_sum(cnt_, c, l);
        return res;
    }
    bool is_sorted(int l, int r) const {  // s[l, r) が広義単調増加（Alphabet の番号の順）か（ABC285 F）  O(log N)
        str_internal::check_range("MutableString::is_sorted", l, r, n_);
        ensure_desc();
        if (r - l <= 1) return true;
        return desc_sum(r - 1) - desc_sum(l) == 0;  // 下がる所 k（s[k] > s[k+1]）が [l, r-1) にない
    }
    long long kth_char(int l, int r, int k) const {  // s[l, r) を並べ替えた k 番目（0 から）の値  O(σ log N)
        str_internal::check_range("MutableString::kth_char", l, r, n_);
        if (k < 0 || k >= r - l) str_internal::fail("MutableString::kth_char", "k は 0〜(r - l - 1)");
        ensure_counts();
        for (int c = 0; c < A_.size(); ++c) {
            const int x = bit_sum(cnt_, c, r) - bit_sum(cnt_, c, l);
            if (k < x) return A_.value(c);
            k -= x;
        }
        return -1;  // ここには来ない
    }
    std::vector<long long> values() const { return val_; }  // 今の値の列
    std::string to_string() const {                         // 今の文字列（値が 0〜255 のとき）
        std::string r(n_, '\0');
        for (int i = 0; i < n_; ++i) {
            if (val_[i] < 0 || val_[i] > 255) str_internal::fail("MutableString::to_string", "値が文字（0〜255）でない");
            r[i] = (char)(unsigned char)val_[i];
        }
        return r;
    }
    const Alphabet& alphabet() const { return A_; }

   private:
    unsigned long long raw(int l, int r) const {  // Σ v_i B^(r-1-i)
        using C = str_internal::RHCore;
        if (l >= r) return 0;
        return C::mul(C::sub(fsum(fw_, r), fsum(fw_, l)), C::pw1[r - 1]);
    }
    unsigned long long raw_rev(int l, int r) const {  // Σ v_i B^(i-l)
        using C = str_internal::RHCore;
        if (l >= r) return 0;
        return C::mul(C::sub(fsum(bw_, r), fsum(bw_, l)), ipw_[l]);
    }
    static unsigned long long fsum(const std::vector<unsigned long long>& t, int i) {  // t の [0, i) の和
        unsigned long long s = 0;
        for (; i > 0; i -= i & -i) s = str_internal::RHCore::add(s, t[i]);
        return s;
    }
    void check_pos(int i, const char* func) const {
        if (i < 0 || i >= n_) str_internal::fail(func, "位置 " + std::to_string(i) + " は 0〜" + std::to_string(n_ - 1) + " の外");
    }
    // 文字ごとのフェニック木（cnt_[c * (N+1) + i]）
    void ensure_counts() const {
        if (counts_built_) return;
        counts_built_ = true;
        const std::size_t W = (std::size_t)n_ + 1;
        if ((long long)A_.size() * (long long)W > 200000000LL) str_internal::fail("MutableString", "σ × (N+1) が 2×10^8（約 800MB）を超える（文字の数の表を作れない）");
        cnt_.assign((std::size_t)A_.size() * W, 0);
        for (int i = 0; i < n_; ++i) ++cnt_[(std::size_t)A_.id(val_[i]) * W + i + 1];
        for (int c = 0; c < A_.size(); ++c) {
            int* t = cnt_.data() + (std::size_t)c * W;
            for (int i = 1; i <= n_; ++i) {
                const int j = i + (i & -i);
                if (j <= n_) t[j] += t[i];
            }
        }
    }
    int bit_sum(const std::vector<int>& t, int c, int i) const {
        if (n_ == 0) return 0;
        const int* b = t.data() + (std::size_t)c * ((std::size_t)n_ + 1);
        int s = 0;
        for (; i > 0; i -= i & -i) s += b[i];
        return s;
    }
    void bit_add(std::vector<int>& t, int c, int i, int d) {
        int* b = t.data() + (std::size_t)c * ((std::size_t)n_ + 1);
        for (int k = i + 1; k <= n_; k += k & -k) b[k] += d;
    }
    // 下がる所（s[k] > s[k+1]）のフェニック木
    void ensure_desc() const {
        if (desc_built_) return;
        desc_built_ = true;
        desc_.assign(n_ + 1, 0), desc_flag_.assign(std::max(0, n_ - 1), 0);
        for (int k = 0; k + 1 < n_; ++k) desc_flag_[k] = A_.id(val_[k]) > A_.id(val_[k + 1]) ? 1 : 0, desc_[k + 1] = desc_flag_[k];
        for (int i = 1; i <= n_; ++i) {
            const int j = i + (i & -i);
            if (j <= n_) desc_[j] += desc_[i];
        }
    }
    int desc_sum(int i) const {
        int s = 0;
        for (; i > 0; i -= i & -i) s += desc_[i];
        return s;
    }
    void bit_add_desc(int k, int d) {
        for (int x = k + 1; x <= n_; x += x & -x) desc_[x] += d;
    }
    Alphabet A_;
    int n_ = 0;
    std::vector<long long> val_;
    std::vector<unsigned long long> ipw_ = std::vector<unsigned long long>(1, 1), fw_ = std::vector<unsigned long long>(1, 0), bw_ = std::vector<unsigned long long>(1, 0);
    mutable std::vector<int> cnt_, desc_;
    mutable std::vector<char> desc_flag_;
    mutable bool counts_built_ = false, desc_built_ = false;
};
#endif  // GG_MUTABLE_STRING
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト
// =============================================================================
namespace ms_test {

#define MS_CHECK(...)                                                                          \
    do {                                                                                        \
        if (!(__VA_ARGS__)) {                                                                   \
            std::cerr << "[mutable_string] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
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

// 愚直：今の列を持ち、毎回その場で計算する
template <class V>
static void check_all(std::mt19937_64& rng, const MutableString& ms, const V& cur, const Alphabet& A) {
    const int n = (int)cur.size();
    MS_CHECK(ms.size() == n);
    std::vector<long long> vals(cur.begin(), cur.end());
    for (auto& x : vals) x = str_internal::code_of((typename V::value_type)x);
    MS_CHECK(ms.values() == vals);
    for (int t = 0; t < 30; ++t) {
        int l = (int)(rng() % (n + 1)), r = (int)(rng() % (n + 1));
        if (l > r) std::swap(l, r);
        const V sub(cur.begin() + l, cur.begin() + r);
        V rev(sub.rbegin(), sub.rend());
        MS_CHECK(ms.hash(l, r) == hash_of(sub) && ms.reversed_hash(l, r) == hash_of(rev));
        MS_CHECK(ms.is_palindrome(l, r) == (sub == rev));
        int l2 = (int)(rng() % (n + 1)), r2 = (int)(rng() % (n + 1));
        if (l2 > r2) std::swap(l2, r2);
        const V sub2(cur.begin() + l2, cur.begin() + r2);
        MS_CHECK(ms.same(l, r, l2, r2) == (sub == sub2));
        MS_CHECK(ms.compare(l, r, l2, r2) == str_internal::lex_compare(sub, sub2));
        int k = 0;
        while (l + k < n && l2 + k < n && cur[l + k] == cur[l2 + k]) ++k;
        MS_CHECK(ms.lcp(l, l2) == k);
        // 文字の数・並び
        std::vector<long long> cnt(A.size(), 0);
        for (auto c : sub) ++cnt[A.id(c)];
        MS_CHECK(ms.counts(l, r) == cnt);
        const int c = (int)(rng() % A.size());
        MS_CHECK(ms.count(l, r, (typename V::value_type)A.value(c)) == cnt[c]);
        bool sorted = true;
        for (int i = l; i + 1 < r; ++i) sorted &= A.id(cur[i]) <= A.id(cur[i + 1]);
        MS_CHECK(ms.is_sorted(l, r) == sorted);
        if (r > l) {
            std::vector<int> ids;
            for (auto x : sub) ids.push_back(A.id(x));
            std::sort(ids.begin(), ids.end());
            const int kk = (int)(rng() % (r - l));
            MS_CHECK(ms.kth_char(l, r, kk) == A.value(ids[kk]));
        }
    }
}

static void test_string(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 300 * TEST_SCALE; ++it) {
        const int sg = 1 + (int)(rng() % 3);
        const Alphabet A = rng() % 2 ? Alphabet::lower() : Alphabet::of("cba");
        std::string s = rnd_str_shaped(rng, (int)(rng() % 30), sg);
        MutableString ms(s, A);
        const bool lazy_first = rng() % 2;  // 表を作る前に書き換える場合も
        for (int q = 0; q < 60; ++q) {
            if (!s.empty() && (q < 10 ? lazy_first : rng() % 2)) {
                const int i = (int)(rng() % s.size());
                const char c = (char)('a' + rng() % sg);
                s[i] = c, ms.set(i, c);
            }
            if (q >= 10 || !lazy_first) check_all(rng, ms, s, A);
        }
        MS_CHECK(ms.to_string() == s);
        ++cases;
    }
}

static void test_values(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 200 * TEST_SCALE; ++it) {
        const Alphabet A = Alphabet::ints(-2, 2);
        std::vector<long long> v(rng() % 25);
        for (auto& x : v) x = (long long)(rng() % 4) - 2;
        MutableString ms(v, A);
        for (int q = 0; q < 40; ++q) {
            if (!v.empty()) {
                const int i = (int)(rng() % v.size());
                v[i] = (long long)(rng() % 4) - 2, ms.set(i, v[i]);
            }
            check_all(rng, ms, v, A);
        }
        ++cases;
    }
    // 大きな値（ハッシュだけ使う。Alphabet の範囲は確かめない）
    for (int it = 0; it < 200 * TEST_SCALE; ++it) {
        std::vector<long long> v(rng() % 25);
        for (auto& x : v) x = rnd_ll(rng);
        MutableString ms(v);
        for (int q = 0; q < 40; ++q) {
            if (!v.empty() && rng() % 2) {
                const int i = (int)(rng() % v.size());
                v[i] = rng() % 2 ? rnd_ll(rng) : v[rng() % v.size()], ms.set(i, v[i]);
            }
            int l = (int)(rng() % (v.size() + 1)), r = (int)(rng() % (v.size() + 1));
            if (l > r) std::swap(l, r);
            const std::vector<long long> sub(v.begin() + l, v.begin() + r);
            MS_CHECK(ms.hash(l, r) == hash_of(sub) && ms.is_palindrome(l, r) == std::equal(sub.begin(), sub.end(), sub.rbegin()));
        }
        ++cases;
    }
}

static void random_test() {
    std::mt19937_64 rng(20261012);
    int cases = 0;
    test_string(rng, cases);
    test_values(rng, cases);
    {  // RollingHash と同じ値
        const std::string s = "abracadabra";
        MutableString ms(s);
        const RollingHash H(s);
        MS_CHECK(ms.hash(2, 9) == H.get(2, 9) && ms.hash(0, 11) == H.all());
        ++cases;
    }
    std::cout << "[mutable_string] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[mutable_string] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    auto ms = [](std::chrono::steady_clock::time_point a) {
        return (long long)std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count();
    };
    const int N = 1000000, Q = 1000000;
    const std::string s = rnd_str(rng, N, 2);
    auto t0 = std::chrono::steady_clock::now();
    MutableString m(s);
    const long long t_build = ms(t0);
    t0 = std::chrono::steady_clock::now();
    long long x = 0;
    for (int q = 0; q < Q; ++q) {
        const int i = (int)(rng() % N);
        if (q % 2) m.set(i, (char)('a' + rng() % 2));
        else x += m.is_palindrome(i, std::min(N, i + (int)(rng() % 20)));
    }
    const long long t_pal = ms(t0);
    t0 = std::chrono::steady_clock::now();
    for (int q = 0; q < 100000; ++q) x += m.lcp((int)(rng() % N), (int)(rng() % N));
    const long long t_lcp = ms(t0);
    t0 = std::chrono::steady_clock::now();
    MutableString m2(rnd_str(rng, 100000, 26));
    for (int q = 0; q < 100000; ++q) {
        const int i = (int)(rng() % 100000);
        if (q % 2) m2.set(i, (char)('a' + rng() % 26));
        else x += m2.is_sorted(i, std::min(100000, i + 3)) + m2.count(0, i, 'c');
    }
    const long long t_cnt = ms(t0);
    std::cout << "[mutable_string] benchmark: 作成（N = 10^6）" << t_build << " ms, set + is_palindrome 10^6 回 " << t_pal << " ms, lcp × 10^5 " << t_lcp
              << " ms, set + is_sorted + count 10^5 回（N = 10^5、表の作成込み）" << t_cnt << " ms  (check " << x % 7 << ")\n";
}

}  // namespace ms_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace ms_demo {

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
#define SHOW(expr, comment) ms_demo::show(#expr, (expr), comment)
#define RUN(...)                   \
    ms_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("書き換えのある文字列（mutable_string） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  set / get / hash / reversed_hash / same / is_palindrome / lcp / compare   書き換え・ハッシュ・比較・回文\n"
              << "  count / counts / is_sorted / kth_char                                    文字の数・並んでいるか・k 番目（表は初めて使うときに作る）\n";

    RUN(MutableString ms(std::string("abcbacb")));

    api("【1】 set(i, c) / get(i) / is_palindrome(l, r)", "なし / long long / bool", "1 文字の書き換えと、区間 [l, r) が回文か（L 文字目〜R 文字目は (L - 1, R)）", "O(log N)");
    {
        note("ABC331 F（Palindrome Query）の入力例 1:");
        std::string out;
        out += ms.is_palindrome(0, 5) ? "Yes " : "No ";
        out += ms.is_palindrome(3, 7) ? "Yes " : "No ";
        out += ms.is_palindrome(1, 2) ? "Yes " : "No ";
        ms.set(4, 'c');
        out += ms.is_palindrome(0, 5) ? "Yes " : "No ";
        out += ms.is_palindrome(3, 7) ? "Yes " : "No ";
        ms.set(3, 'c');
        out += ms.is_palindrome(2, 6) ? "Yes " : "No ";
        note("  → " + out);
        SHOW(ms.to_string(), "");
        SHOW((char)ms.get(4), "");
    }

    api("【2】 hash(l, r) / reversed_hash(l, r) / same(l1, r1, l2, r2)", "StrHash / StrHash / bool",
        "区間のハッシュ（01_rolling_hash の RollingHash::get・hash_of と同じ値）/ 反転したもののハッシュ / 区間どうしが等しいか", "O(log N)");
    SHOW(ms.same(1, 3, 4, 6), "\"bc\" と \"cc\"");
    SHOW(ms.hash(0, 3) == hash_of("abc"), "ほかの文字列のハッシュとも比べられる");
    SHOW(ms.reversed_hash(0, 3) == hash_of("cba"), "");

    api("【3】 lcp(i, j) / compare(l1, r1, l2, r2)", "int / int（-1 / 0 / 1）", "書き換えたあとの最長共通接頭辞・辞書順（値の大小で）", "O(log^2 N)");
    SHOW(ms.lcp(2, 3), "");
    SHOW(ms.compare(0, 2, 5, 7), "\"ab\" と \"cb\"");

    api("【4】 count(l, r, c) / counts(l, r) / is_sorted(l, r) / kth_char(l, r, k)", "int / vector<int> / bool / long long",
        "区間の文字 c の数・番号ごとの数・広義単調増加か（Alphabet の番号の順）・並べ替えた k 番目（0 から）", "O(log N)・O(σ log N)");
    {
        note("ABC285 F（Substring of Sorted String）: s[l, r) が並べ替えた S の部分文字列 ⇔ is_sorted(l, r) かつ");
        note("  両端の文字より真ん中の文字（s[l] < c < s[r-1]）がすべて区間の中にある（count(l, r, c) == count(0, N, c)）");
        MutableString t(std::string("abcdcf"));
        auto ok = [&](int l, int r) {
            if (!t.is_sorted(l, r)) return false;
            const int a = t.alphabet().id(t.get(l)), b = t.alphabet().id(t.get(r - 1));
            for (int c = a + 1; c < b; ++c)
                if (t.count(l, r, t.alphabet().ch(c)) != t.count(0, t.size(), t.alphabet().ch(c))) return false;
            return true;
        };
        std::string out;
        out += ok(0, 3) ? "Yes " : "No ";
        out += ok(1, 6) ? "Yes " : "No ";
        t.set(4, 'e');
        out += ok(1, 6) ? "Yes " : "No ";
        note("  入力例 1 → " + out);
        SHOW(t.counts(0, 6), "a〜z の数");
        SHOW((char)t.kth_char(0, 6, 3), "abcdef の 3 番目");
    }
}

}  // namespace ms_demo

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
    ms_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    ms_test::benchmark();    // 速度計測
    ms_demo::demo();         // 使い方デモ
    return 0;
}
