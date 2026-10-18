// =============================================================================
//  Aho–Corasick（多数のパターンの照合：出現回数・最初の出現・すべての出現・含むか・遷移の表（DP 用）・パターンどうしの包含）   --- 文字列第二世代
// -----------------------------------------------------------------------------
//  分類       : 02 一致・出現（パターン検索）
//  依存       : 00_str_base（同梱済み）
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）
//  文字コード : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  AhoCorasick ac(A)：add(p) でパターンを足し（番号 0, 1, … を返す。同じ文字列を 2 回足してもよい）、build() してから使う
//  状態 = Trie のノード（0 が根）。遷移の表は ノード数 × σ 個の int（パターンの長さの合計 10^6・σ = 26 で約 100MB）
//  パターン 1 つなら 02_pattern_match（find_all・KmpAutomaton）で足りる
//
//  機能一覧                                                                       計算量
//    ac.add(p) / build()                   パターンを足す（番号を返す）/ 遷移と失敗のリンクを作る           O(|p|) / O(ノード数 × σ)
//    ac.count_each(t)                      パターンごとの出現回数（ABC362 G）                             O(|t| + ノード数)
//    ac.count_total(t) / contains_any(t)   出現回数の合計 / どれかが現れるか                               O(|t|)
//    ac.first_occurrence_each(t)           パターンごとの最初の出現の開始位置（なければ -1）                 O(|t| + ノード数)
//    ac.for_each_match(t, f)               すべての出現で f(開始位置, 番号)                                O(|t| + 出現の数)
//    ac.contained_patterns(id)             パターン id の中に現れるほかのパターン（包含関係。ABC354 G）       O(|p| + 個数)
//    状態（DP 用）: next(s, c) / next_id(s, id)（O(1)）/ match_count(s) / longest_match(s) / match_ids(s)
//                  depth(s) / link(s) / node_of(id) / pattern_length(id) / bfs_order() / size() / pattern_count()
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    ac_test::random_test()  … すべての位置でパターンを比べる愚直解と全機能を突き合わせるランダムテスト
//    ac_test::benchmark()    … 速度計測
//    ac_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
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
#ifndef GG_AHO_CORASICK
#define GG_AHO_CORASICK
// ---- Aho–Corasick（多数のパターンの照合。文字列第二世代） ----
// AhoCorasick ac(A)：add(p) でパターンを足し（番号 0, 1, … を返す。同じ文字列を 2 回足してもよい）、build() してから使う
// 状態 = Trie のノード（0 が根）。next(state, c) は O(1)（遷移の表は ノード数 × σ 個の int。2×10^8 を超えると止める）。build の後に add はできない
// 本文・パターンの文字は Alphabet の範囲で（範囲外なら止める）
struct AhoCorasick {
    explicit AhoCorasick(const Alphabet& A = Alphabet::lower()) : A_(A), sg_(A.size()) { new_node(-1, 0); }

    template <class Seq>
    int add(const Seq& p0) {  // パターン p を足して番号を返す（空のパターンは止める）  O(|p|)
        if (built_) str_internal::fail("AhoCorasick::add", "build() の後には足せない");
        const auto& p = str_internal::seq(p0);
        if (p.size() == 0) str_internal::fail("AhoCorasick::add", "空のパターンは使えない");
        int v = 0;
        for (std::size_t i = 0; i < p.size(); ++i) {
            const int c = A_.id(p[i]);
            int u = go_[(std::size_t)v * sg_ + c];
            if (u < 0) u = new_node(v, dep_[v] + 1), go_[(std::size_t)v * sg_ + c] = u;
            v = u;
        }
        node_of_.push_back(v), len_.push_back((int)p.size());
        return (int)node_of_.size() - 1;
    }
    void build() {  // 失敗のリンクと遷移の表を作る  O(ノード数 × σ)
        if (built_) return;
        built_ = true;
        const int n = size();
        link_.assign(n, 0), out_.assign(n, -1), cnt_.assign(n, 0), longest_.assign(n, 0), order_.clear();
        // ノードごとのパターンの番号（CSR）
        ids_start_.assign(n + 1, 0);
        for (int v : node_of_) ++ids_start_[v + 1], ++cnt_[v];
        for (int v = 0; v < n; ++v) ids_start_[v + 1] += ids_start_[v];
        ids_.assign(node_of_.size(), 0);
        std::vector<int> fill(ids_start_.begin(), ids_start_.end() - 1);
        for (int id = 0; id < (int)node_of_.size(); ++id) ids_[fill[node_of_[id]]++] = id;
        // 幅優先
        order_.reserve(n);
        order_.push_back(0);
        for (int c = 0; c < sg_; ++c) {
            int& u = go_[c];
            if (u < 0) u = 0;
            else link_[u] = 0, order_.push_back(u);
        }
        for (std::size_t h = 1; h < order_.size(); ++h) {
            const int v = order_[h];
            for (int c = 0; c < sg_; ++c) {
                int& u = go_[(std::size_t)v * sg_ + c];
                const int f = go_[(std::size_t)link_[v] * sg_ + c];
                if (u < 0) u = f;
                else link_[u] = f, order_.push_back(u);
            }
        }
        for (int v : order_) {  // 根から順に：リンクの先までの数・最長・出力のリンク
            if (v == 0) continue;
            const int f = link_[v];
            out_[v] = (ids_start_[f + 1] > ids_start_[f]) ? f : out_[f];
            longest_[v] = (ids_start_[v + 1] > ids_start_[v]) ? dep_[v] : longest_[f];
            cnt_[v] += cnt_[f];
        }
    }

    // ---- 状態（DP 用） ----
    int size() const { return (int)dep_.size(); }                // 状態の数
    int pattern_count() const { return (int)node_of_.size(); }   // パターンの数
    template <class C>
    int next(int state, C c) const {  // 状態 state で文字 c を読んだ次の状態  O(1)
        check_built("AhoCorasick::next"), check_state(state, "AhoCorasick::next");
        return go_[(std::size_t)state * sg_ + A_.id(c)];
    }
    int next_id(int state, int id) const {  // 文字の番号 id で  O(1)
        check_built("AhoCorasick::next_id"), check_state(state, "AhoCorasick::next_id");
        if (id < 0 || id >= sg_) str_internal::fail("AhoCorasick::next_id", "番号は 0〜σ-1");
        return go_[(std::size_t)state * sg_ + id];
    }
    long long match_count(int state) const { return check_built("AhoCorasick::match_count"), check_state(state, "AhoCorasick::match_count"), cnt_[state]; }  // ここで終わるパターンの数（同じ文字列も別に数える）
    int longest_match(int state) const { return check_built("AhoCorasick::longest_match"), check_state(state, "AhoCorasick::longest_match"), longest_[state]; }  // ここで終わる最長のパターンの長さ（なければ 0）
    std::vector<int> match_ids(int state) const {  // ここで終わるパターンの番号（長い順）  O(個数)
        check_built("AhoCorasick::match_ids"), check_state(state, "AhoCorasick::match_ids");
        std::vector<int> r;
        for (int v = (ids_start_[state + 1] > ids_start_[state]) ? state : out_[state]; v > 0; v = out_[v])
            for (int k = ids_start_[v]; k < ids_start_[v + 1]; ++k) r.push_back(ids_[k]);
        return r;
    }
    int depth(int state) const { return check_state(state, "AhoCorasick::depth"), dep_[state]; }                                   // 根からの深さ（状態の文字列の長さ）
    int link(int state) const { return check_built("AhoCorasick::link"), check_state(state, "AhoCorasick::link"), link_[state]; }  // 失敗のリンク（根は 0）
    int node_of(int id) const { return check_id(id, "AhoCorasick::node_of"), node_of_[id]; }  // パターン id が終わる状態
    int pattern_length(int id) const { return check_id(id, "AhoCorasick::pattern_length"), len_[id]; }
    const std::vector<int>& bfs_order() const { return check_built("AhoCorasick::bfs_order"), order_; }  // 根から幅優先の順（リンクの先が先に来る）

    // ---- 文字列の中の出現 ----
    template <class Seq>
    std::vector<long long> count_each(const Seq& t0) const {  // パターンごとの出現回数（ABC362 G）  O(|t| + ノード数)
        check_built("AhoCorasick::count_each");
        const auto& t = str_internal::seq(t0);
        std::vector<long long> vis(size(), 0);
        int v = 0;
        for (std::size_t i = 0; i < t.size(); ++i) v = go_[(std::size_t)v * sg_ + A_.id(t[i])], ++vis[v];
        for (int k = size() - 1; k > 0; --k) vis[link_[order_[k]]] += vis[order_[k]];
        std::vector<long long> r(pattern_count());
        for (int id = 0; id < pattern_count(); ++id) r[id] = vis[node_of_[id]];
        return r;
    }
    template <class Seq>
    long long count_total(const Seq& t0) const {  // すべてのパターンの出現回数の合計  O(|t|)
        check_built("AhoCorasick::count_total");
        const auto& t = str_internal::seq(t0);
        long long s = 0;
        int v = 0;
        for (std::size_t i = 0; i < t.size(); ++i) v = go_[(std::size_t)v * sg_ + A_.id(t[i])], s += cnt_[v];
        return s;
    }
    template <class Seq>
    std::vector<int> first_occurrence_each(const Seq& t0) const {  // パターンごとの最初の出現の開始位置（なければ -1）  O(|t| + ノード数)
        check_built("AhoCorasick::first_occurrence_each");
        const auto& t = str_internal::seq(t0);
        const int INFPOS = std::numeric_limits<int>::max();
        std::vector<int> first(size(), INFPOS);
        int v = 0;
        for (std::size_t i = 0; i < t.size(); ++i) {
            v = go_[(std::size_t)v * sg_ + A_.id(t[i])];
            if (first[v] == INFPOS) first[v] = (int)i;
        }
        for (int k = size() - 1; k > 0; --k) first[link_[order_[k]]] = std::min(first[link_[order_[k]]], first[order_[k]]);
        std::vector<int> r(pattern_count());
        for (int id = 0; id < pattern_count(); ++id) r[id] = first[node_of_[id]] == INFPOS ? -1 : first[node_of_[id]] - len_[id] + 1;
        return r;
    }
    template <class Seq, class F>
    void for_each_match(const Seq& t0, F f) const {  // すべての出現で f(開始位置, パターンの番号)（終わる位置の順）  O(|t| + 出現の数)
        check_built("AhoCorasick::for_each_match");
        const auto& t = str_internal::seq(t0);
        int v = 0;
        for (std::size_t i = 0; i < t.size(); ++i) {
            v = go_[(std::size_t)v * sg_ + A_.id(t[i])];
            for (int u = (ids_start_[v + 1] > ids_start_[v]) ? v : out_[v]; u > 0; u = out_[u])
                for (int k = ids_start_[u]; k < ids_start_[u + 1]; ++k) f((int)i - len_[ids_[k]] + 1, ids_[k]);
        }
    }
    template <class Seq>
    bool contains_any(const Seq& t0) const {  // どれかのパターンが現れるか  O(|t|)
        check_built("AhoCorasick::contains_any");
        const auto& t = str_internal::seq(t0);
        int v = 0;
        for (std::size_t i = 0; i < t.size(); ++i) {
            v = go_[(std::size_t)v * sg_ + A_.id(t[i])];
            if (cnt_[v] > 0) return true;
        }
        return false;
    }
    std::vector<int> contained_patterns(int id) const {  // パターン id の中に（部分文字列として）現れるほかのパターンの番号（昇順。同じ文字列のパターンも入る）  O(|p| + 個数 log 個数)
        check_built("AhoCorasick::contained_patterns");
        check_id(id, "AhoCorasick::contained_patterns");
        // id の接頭辞の状態それぞれから出力のリンクをたどって、そこで終わるパターンを集める（一度見た状態から先は集めたあとなので止める）
        if (mark_.size() != (std::size_t)size()) mark_.assign(size(), 0);
        ++stamp_;
        std::vector<int> r;
        for (int v : path_nodes(id)) {
            for (int u = (ids_start_[v + 1] > ids_start_[v]) ? v : out_[v]; u > 0 && mark_[u] != stamp_; u = out_[u]) {
                mark_[u] = stamp_;
                for (int k = ids_start_[u]; k < ids_start_[u + 1]; ++k)
                    if (ids_[k] != id) r.push_back(ids_[k]);
            }
        }
        std::sort(r.begin(), r.end());
        return r;
    }
    const Alphabet& alphabet() const { return A_; }

   private:
    int new_node(int parent, int d) {
        if ((long long)(dep_.size() + 1) * sg_ > 200000000LL) str_internal::fail("AhoCorasick", "ノード数 × σ が 2×10^8（遷移の表が約 800MB）を超える（σ が大きい数の列には向かない）");
        go_.insert(go_.end(), (std::size_t)sg_, -1);
        dep_.push_back(d), par_.push_back(parent);
        return (int)dep_.size() - 1;
    }
    std::vector<int> path_nodes(int id) const {  // パターン id の接頭辞の状態（長さ 1〜|p|）
        std::vector<int> r(len_[id]);
        int v = node_of_[id];
        for (int k = len_[id] - 1; k >= 0; --k) r[k] = v, v = par_[v];
        return r;
    }
    void check_built(const char* func) const {
        if (!built_) str_internal::fail(func, "先に build() を呼ぶ");
    }
    void check_state(int s, const char* func) const {
        if (s < 0 || s >= size()) str_internal::fail(func, "状態 " + std::to_string(s) + " は 0〜size()-1 の外");
    }
    void check_id(int id, const char* func) const {
        if (id < 0 || id >= pattern_count()) str_internal::fail(func, "パターンの番号 " + std::to_string(id) + " は 0〜pattern_count()-1 の外");
    }
    Alphabet A_;
    int sg_;
    bool built_ = false;
    std::vector<int> go_, dep_, par_, node_of_, len_, link_, out_, longest_, order_, ids_start_, ids_;
    std::vector<long long> cnt_;
    mutable std::vector<int> mark_;  // contained_patterns の「見た」印
    mutable int stamp_ = 0;
};
#endif  // GG_AHO_CORASICK
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト
// =============================================================================
namespace ac_test {

#define AC_CHECK(...)                                                                          \
    do {                                                                                        \
        if (!(__VA_ARGS__)) {                                                                   \
            std::cerr << "[aho_corasick] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
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

static bool ends_with_at(const std::string& t, int end, const std::string& p) {  // t[0, end) の末尾が p か
    return (int)p.size() <= end && t.compare(end - p.size(), p.size(), p) == 0;
}

static void test_basic(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 1500 * TEST_SCALE; ++it) {
        const int sg = 1 + (int)(rng() % 3);
        const int k = 1 + (int)(rng() % 6);
        std::vector<std::string> ps(k);
        for (auto& p : ps) p = rnd_str_shaped(rng, 1 + (int)(rng() % 5), sg);
        if (k >= 2 && rng() % 4 == 0) ps[1] = ps[0];  // 同じ文字列
        const std::string t = rnd_str_shaped(rng, (int)(rng() % 40), sg);
        AhoCorasick ac;
        for (int i = 0; i < k; ++i) AC_CHECK(ac.add(ps[i]) == i);
        ac.build();
        AC_CHECK(ac.pattern_count() == k);
        // 出現
        std::vector<long long> cnt(k, 0);
        std::vector<int> first(k, -1);
        std::vector<std::pair<int, int>> occ;  // (終わり, 番号) の順に並べる
        long long total = 0;
        for (int e = 1; e <= (int)t.size(); ++e)
            for (int i = 0; i < k; ++i)
                if (ends_with_at(t, e, ps[i])) {
                    ++cnt[i], ++total;
                    if (first[i] < 0) first[i] = e - (int)ps[i].size();
                    occ.push_back({e, i});
                }
        AC_CHECK(ac.count_each(t) == cnt && ac.count_total(t) == total && ac.first_occurrence_each(t) == first && ac.contains_any(t) == (total > 0));
        std::vector<std::pair<int, int>> got;
        ac.for_each_match(t, [&](int st, int id) { got.push_back({st + (int)ps[id].size(), id}); });
        std::vector<std::pair<int, int>> a = got, b = occ;
        std::sort(a.begin(), a.end()), std::sort(b.begin(), b.end());
        AC_CHECK(a == b);
        for (std::size_t q = 1; q < got.size(); ++q) AC_CHECK(got[q - 1].first <= got[q].first);  // 終わる位置の順
        // 状態：t の各接頭辞を読んだあと
        int s = 0;
        for (int e = 1; e <= (int)t.size(); ++e) {
            const int s2 = ac.next(s, t[e - 1]);
            AC_CHECK(s2 == ac.next_id(s, ac.alphabet().id(t[e - 1])));
            s = s2;
            long long mc = 0;
            int lm = 0, dep = 0;
            std::vector<int> ids;
            for (int i = 0; i < k; ++i)
                if (ends_with_at(t, e, ps[i])) ++mc, lm = std::max(lm, (int)ps[i].size()), ids.push_back(i);
            for (int d = 0; d <= e; ++d) {  // t[0, e) の接尾辞で、どれかのパターンの接頭辞であるものの最長
                const std::string suf = t.substr(e - d, d);
                for (const auto& p : ps)
                    if (p.compare(0, d, suf) == 0 && (int)p.size() >= d) dep = std::max(dep, d);
            }
            AC_CHECK(ac.match_count(s) == mc && ac.longest_match(s) == lm && ac.depth(s) == dep);
            auto mi = ac.match_ids(s);
            std::sort(mi.begin(), mi.end());
            AC_CHECK(mi == ids);
        }
        // 包含
        for (int i = 0; i < k; ++i) {
            std::vector<int> ex;
            for (int j = 0; j < k; ++j)
                if (j != i && ps[i].find(ps[j]) != std::string::npos) ex.push_back(j);
            AC_CHECK(ac.contained_patterns(i) == ex);
            AC_CHECK(ac.depth(ac.node_of(i)) == (int)ps[i].size() && ac.pattern_length(i) == (int)ps[i].size());
        }
        // 幅優先の順でリンクの先が先
        std::vector<int> pos(ac.size(), -1);
        const auto& ord = ac.bfs_order();
        AC_CHECK((int)ord.size() == ac.size() && ord[0] == 0);
        for (int q = 0; q < (int)ord.size(); ++q) pos[ord[q]] = q;
        for (int v = 1; v < ac.size(); ++v) AC_CHECK(pos[ac.link(v)] < pos[v] && ac.depth(ac.link(v)) < ac.depth(v));
        ++cases;
    }
}

static void test_ints(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 500 * TEST_SCALE; ++it) {
        const Alphabet A = Alphabet::ints(-1, 2);
        const int k = 1 + (int)(rng() % 4);
        std::vector<std::vector<int>> ps(k);
        for (auto& p : ps) {
            p.resize(1 + rng() % 4);
            for (auto& x : p) x = (int)(rng() % 3) - 1;
        }
        std::vector<long long> t(rng() % 30);
        for (auto& x : t) x = (long long)(rng() % 3) - 1;
        AhoCorasick ac(A);
        for (auto& p : ps) ac.add(p);
        ac.build();
        std::vector<long long> cnt(k, 0);
        for (int i = 0; i < k; ++i)
            for (int st = 0; st + (int)ps[i].size() <= (int)t.size(); ++st) cnt[i] += std::equal(ps[i].begin(), ps[i].end(), t.begin() + st);
        AC_CHECK(ac.count_each(t) == cnt);
        ++cases;
    }
}

static void random_test() {
    std::mt19937_64 rng(20261012);
    int cases = 0;
    test_basic(rng, cases);
    test_ints(rng, cases);
    std::cout << "[aho_corasick] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[aho_corasick] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    auto ms = [](std::chrono::steady_clock::time_point a) {
        return (long long)std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count();
    };
    const std::string t = rnd_str(rng, 1000000, 4);
    std::vector<std::string> ps;
    long long tot = 0;
    while (tot < 1000000) {
        const int len = 1 + (int)(rng() % 12), st = (int)(rng() % (t.size() - len));
        ps.push_back(t.substr(st, len)), tot += len;
    }
    auto t0 = std::chrono::steady_clock::now();
    AhoCorasick ac;
    for (const auto& p : ps) ac.add(p);
    ac.build();
    const long long t_build = ms(t0);
    t0 = std::chrono::steady_clock::now();
    const auto c = ac.count_each(t);
    const long long t_cnt = ms(t0);
    long long x = c[0] + ac.count_total(t);
    t0 = std::chrono::steady_clock::now();
    AhoCorasick big;
    std::vector<std::string> ps2;
    tot = 0;
    while (tot < 1000000) ps2.push_back(rnd_str(rng, 1 + (int)(rng() % 20), 26)), tot += (long long)ps2.back().size();
    for (const auto& p : ps2) big.add(p);
    big.build();
    x += big.count_each(rnd_str(rng, 1000000, 26))[0];
    const long long t_big = ms(t0);
    std::cout << "[aho_corasick] benchmark: パターンの長さの合計 10^6（σ = 4）build " << t_build << " ms, count_each（|t| = 10^6）" << t_cnt
              << " ms, σ = 26 で build + count_each " << t_big << " ms, ノード " << big.size() << "  (check " << x % 7 << ")\n";
}

}  // namespace ac_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace ac_demo {

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
#define SHOW(expr, comment) ac_demo::show(#expr, (expr), comment)
#define RUN(...)                   \
    ac_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("Aho–Corasick（aho_corasick） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  add / build                                          パターンを足す・作る\n"
              << "  count_each / count_total / contains_any / first_occurrence_each / for_each_match   文字列の中の出現\n"
              << "  contained_patterns                                   パターンどうしの包含\n"
              << "  next / next_id / match_count / longest_match / match_ids / depth / link / bfs_order   状態（DP 用）\n";

    RUN(AhoCorasick ac);
    RUN(for (const std::string p : {"i", "s", "a", "is", "missisippi"}) ac.add(p));
    RUN(ac.build());

    api("【1】 count_each(t) / count_total(t) / contains_any(t)", "vector<long long> / long long / bool", "パターンごとの出現回数 / その合計 / どれかが現れるか", "O(|t| + ノード数)");
    SHOW(ac.count_each("missisippi"), "ABC362 G（Count Substring Query）の入力例 1");
    SHOW(ac.count_total("missisippi"), "");
    SHOW(ac.contains_any("xyz"), "");

    api("【2】 first_occurrence_each(t) / for_each_match(t, f)", "vector<int> / なし（f(開始位置, 番号)）", "パターンごとの最初の出現（なければ -1）/ すべての出現（終わる位置の順）", "O(|t| + ノード数) / O(|t| + 出現の数)");
    SHOW(ac.first_occurrence_each("misa"), "");
    {
        std::string out;
        ac.for_each_match("misis", [&](int st, int id) { out += "(" + std::to_string(st) + "," + std::to_string(id) + ")"; });
        code("ac.for_each_match(\"misis\", f) で (開始位置, 番号) -> " + out);
    }

    api("【3】 contained_patterns(id)", "vector<int>", "パターン id の中に部分文字列として現れるほかのパターン（同じ文字列も）。包含関係の DAG・反鎖（ABC354 G は最大流で）", "O(|p| + 個数)");
    {
        AhoCorasick a2;
        for (const std::string p : {"atcoder", "at", "coder", "code"}) a2.add(p);
        a2.build();
        code("{\"atcoder\", \"at\", \"coder\", \"code\"}   // ABC354 G（Select Strings）の入力例 1");
        SHOW(a2.contained_patterns(0), "atcoder の中に at・coder・code");
        SHOW(a2.contained_patterns(2), "coder の中に code");
        note("同じ文字列は互いに含むので、ABC354 G では先に重複を除く（重みの大きい方を残す）");
    }

    api("【4】 状態：next(s, c) / match_count(s) / longest_match(s) / match_ids(s) / depth(s)", "int / long long / int / vector<int> / int",
        "読んだ文字列の「どれかのパターンの接頭辞である最長の接尾辞」が状態。そこで終わるパターンの数・最長・番号", "O(1)（match_ids は個数）");
    {
        int s = 0;
        for (char c : std::string("mis")) s = ac.next(s, c);
        code("\"mis\" を読んだ状態 s");
        SHOW(ac.match_count(s), "\"mis\" の末尾で終わる：is と s");
        SHOW(ac.depth(s), "\"mis\" は missisippi の接頭辞");
        s = ac.next(s, 'i'), s = ac.next(s, 's');
        SHOW(ac.match_ids(s), "\"misis\" の末尾で終わる：is, s（長い順）");
        SHOW(ac.longest_match(s), "");
    }
    {
        note("ABC268 Ex（Taboo）: 状態を進め、match_count(s) > 0 になったらその文字を * にして根（0）に戻す（貪欲）");
        AhoCorasick a3;
        for (const std::string p : {"abcd", "ijk", "ghi"}) a3.add(p);
        a3.build();
        int s = 0, ans = 0;
        for (char c : std::string("abcdefghijklmn")) {
            s = a3.next(s, c);
            if (a3.match_count(s) > 0) ++ans, s = 0;
        }
        note("  入力例 1 → " + to_s(ans));
    }
    {
        note("パターンを含まない長さ L の文字列の数（DP）: dp[状態] を next で進め、match_count(s) > 0 の状態へは行かない");
        AhoCorasick a4(Alphabet::of("ab"));
        a4.add("aa"), a4.add("bab");
        a4.build();
        std::vector<long long> dp(a4.size(), 0);
        dp[0] = 1;
        for (int step = 0; step < 6; ++step) {
            std::vector<long long> nd(a4.size(), 0);
            for (int v = 0; v < a4.size(); ++v)
                for (int c = 0; c < 2; ++c) {
                    const int u = a4.next_id(v, c);
                    if (a4.match_count(u) == 0) nd[u] += dp[v];
                }
            dp = nd;
        }
        long long total = 0;
        for (long long x : dp) total += x;
        note("  {aa, bab} を含まない a・b の長さ 6 の文字列 → " + to_s(total) + " 個（L が 10^18 なら math_v2 05_matrix で遷移行列を累乗）");
    }
}

}  // namespace ac_demo

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
    ac_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    ac_test::benchmark();    // 速度計測
    ac_demo::demo();         // 使い方デモ
    return 0;
}
