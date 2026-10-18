// =============================================================================
//  ワイルドカードつきの一致・ずらしごとの不一致の数・k 個以下の不一致で一致（畳み込み NTT）   --- 文字列第二世代
// -----------------------------------------------------------------------------
//  分類       : 02 一致・出現（パターン検索）
//  依存       : 00_str_base（同梱済み）
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）
//  文字コード : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  中で NTT（mod 998244353）を使う（math_v2 には依存しない）。text の長さ N、pat の長さ M。ずらし i は 0〜N-M（pat を text[i, i+M) に重ねる）
//  wild は text・pat のどちらにあっても任意の 1 文字と一致する。列は string でも数の列でもよい（数の列で wild を使うときは値を渡す）
//  wild のない普通の一致は 02_pattern_match の find_all（O(N + M) で速い）
//
//  機能一覧                                                                       計算量
//    mismatch_counts(text, pat, wild)      ずらしごとの不一致の数（wild を省くとワイルドカードなし）         O(σ' (N + M) log(N + M))
//    match_with_k_mismatches(text, pat, k, wild)   不一致が k 以下のずらしの一覧                         同上
//    wildcard_match(text, pat, wild = '?') wild を任意の 1 文字として一致するずらしの一覧（乱数を使う）      O((N + M) log(N + M))
//    （σ' は pat に出てくる異なる値の数。値の種類が多いと遅い）
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    wc_test::random_test()  … すべてのずらしで 1 文字ずつ比べる愚直解と全機能を突き合わせるランダムテスト
//    wc_test::benchmark()    … 速度計測
//    wc_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
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
#ifndef GG_WILDCARD_MATCH
#define GG_WILDCARD_MATCH
// ---- ワイルドカードつきの一致・ずらしごとの不一致の数（畳み込み。文字列第二世代） ----
// 中で NTT（mod 998244353）を使う（math_v2 には依存しない）。text の長さ N（2^23 まで）、pat の長さ M。ずらし i は 0〜N-M
// wild（文字列の wildcard_match では省くと '?'）は text・pat のどちらにあっても任意の 1 文字と一致する
namespace str_internal {
struct WcNtt {  // mod 998244353 の NTT（4 基数。結果はビットを逆にした順で、点ごとの積と ntt_inv の組で使う）
    static constexpr unsigned MOD = 998244353;
    static unsigned mul(unsigned a, unsigned b) { return (unsigned)((unsigned long long)a * b % MOD); }
    static unsigned add(unsigned a, unsigned b) { return a + b >= MOD ? a + b - MOD : a + b; }
    static unsigned sub(unsigned a, unsigned b) { return a >= b ? a - b : a + MOD - b; }
    static int size_for(long long need, const char* func) {  // need 以上の最小の 2 冪（2^23 まで）
        if (need > (1LL << 23)) fail(func, "text の長さ N が 2^23 を超える（NTT の長さの上限）");
        int sz = 1;
        while (sz < need) sz <<= 1;
        return sz;
    }
    static unsigned pw(unsigned long long a, unsigned long long e) {
        unsigned long long r = 1;
        a %= MOD;
        for (; e; e >>= 1, a = a * a % MOD)
            if (e & 1) r = r * a % MOD;
        return (unsigned)r;
    }
    struct Table {
        unsigned root[24], iroot[24], rate2[23], irate2[23], rate3[22], irate3[22];
        Table() {
            const int r = 23;  // 998244353 - 1 = 119 × 2^23
            root[r] = pw(3, (MOD - 1) >> r), iroot[r] = pw(root[r], MOD - 2);
            for (int i = r - 1; i >= 0; --i) root[i] = mul(root[i + 1], root[i + 1]), iroot[i] = mul(iroot[i + 1], iroot[i + 1]);
            unsigned prod = 1, iprod = 1;
            for (int i = 0; i + 2 <= r; ++i) {
                rate2[i] = mul(root[i + 2], prod), irate2[i] = mul(iroot[i + 2], iprod);
                prod = mul(prod, iroot[i + 2]), iprod = mul(iprod, root[i + 2]);
            }
            prod = 1, iprod = 1;
            for (int i = 0; i + 3 <= r; ++i) {
                rate3[i] = mul(root[i + 3], prod), irate3[i] = mul(iroot[i + 3], iprod);
                prod = mul(prod, iroot[i + 3]), iprod = mul(iprod, root[i + 3]);
            }
        }
    };
    static const Table& table() {
        static const Table t;
        return t;
    }
    static void ntt(std::vector<unsigned>& a) {  // 長さ 2^h
        const int n = (int)a.size();
        int h = 0;
        while ((1 << h) < n) ++h;
        const Table& t = table();
        for (int len = 0; len < h;) {
            if (h - len == 1) {
                const int p = 1 << (h - len - 1);
                unsigned rot = 1;
                for (int s = 0; s < (1 << len); ++s) {
                    const int off = s << (h - len);
                    for (int i = 0; i < p; ++i) {
                        const unsigned l = a[i + off], r = mul(a[i + off + p], rot);
                        a[i + off] = add(l, r), a[i + off + p] = sub(l, r);
                    }
                    if (s + 1 != (1 << len)) rot = mul(rot, t.rate2[__builtin_ctz(~(unsigned)s)]);
                }
                ++len;
            } else {
                const int p = 1 << (h - len - 2);
                unsigned rot = 1;
                const unsigned imag = t.root[2];
                for (int s = 0; s < (1 << len); ++s) {
                    const unsigned rot2 = mul(rot, rot), rot3 = mul(rot2, rot);
                    const int off = s << (h - len);
                    for (int i = 0; i < p; ++i) {
                        const unsigned a0 = a[i + off], a1 = mul(a[i + off + p], rot), a2 = mul(a[i + off + 2 * p], rot2), a3 = mul(a[i + off + 3 * p], rot3);
                        const unsigned x = mul(sub(a1, a3), imag), s02 = add(a0, a2), d02 = sub(a0, a2), s13 = add(a1, a3);
                        a[i + off] = add(s02, s13);
                        a[i + off + p] = sub(s02, s13);
                        a[i + off + 2 * p] = add(d02, x);
                        a[i + off + 3 * p] = sub(d02, x);
                    }
                    if (s + 1 != (1 << len)) rot = mul(rot, t.rate3[__builtin_ctz(~(unsigned)s)]);
                }
                len += 2;
            }
        }
    }
    static void ntt_inv(std::vector<unsigned>& a) {  // ntt の逆（1/n 倍まで）
        const int n = (int)a.size();
        int h = 0;
        while ((1 << h) < n) ++h;
        const Table& t = table();
        for (int len = h; len;) {
            if (len == 1) {
                const int p = 1 << (h - len);
                unsigned irot = 1;
                for (int s = 0; s < (1 << (len - 1)); ++s) {
                    const int off = s << (h - len + 1);
                    for (int i = 0; i < p; ++i) {
                        const unsigned l = a[i + off], r = a[i + off + p];
                        a[i + off] = add(l, r), a[i + off + p] = mul(sub(l, r), irot);
                    }
                    if (s + 1 != (1 << (len - 1))) irot = mul(irot, t.irate2[__builtin_ctz(~(unsigned)s)]);
                }
                --len;
            } else {
                const int p = 1 << (h - len);
                unsigned irot = 1;
                const unsigned iimag = t.iroot[2];
                for (int s = 0; s < (1 << (len - 2)); ++s) {
                    const unsigned irot2 = mul(irot, irot), irot3 = mul(irot2, irot);
                    const int off = s << (h - len + 2);
                    for (int i = 0; i < p; ++i) {
                        const unsigned a0 = a[i + off], a1 = a[i + off + p], a2 = a[i + off + 2 * p], a3 = a[i + off + 3 * p];
                        const unsigned x = mul(sub(a2, a3), iimag);
                        a[i + off] = add(add(a0, a1), add(a2, a3));
                        a[i + off + p] = mul(add(sub(a0, a1), x), irot);
                        a[i + off + 2 * p] = mul(sub(add(a0, a1), add(a2, a3)), irot2);
                        a[i + off + 3 * p] = mul(sub(sub(a0, a1), x), irot3);
                    }
                    if (s + 1 != (1 << (len - 2))) irot = mul(irot, t.irate3[__builtin_ctz(~(unsigned)s)]);
                }
                len -= 2;
            }
        }
        const unsigned inv_n = pw(n, MOD - 2);
        for (auto& x : a) x = mul(x, inv_n);
    }
};
// 本体（has_wild = false なら wild は使わない）
template <class T, class P>
std::vector<int> mismatch_counts_impl(const T& text0, const P& pat0, bool has_wild, long long wild) {
    using W = WcNtt;
    const auto t = codes_of(text0), p = codes_of(pat0);
    const int n = (int)t.size(), m = (int)p.size();
    if (m > n) return {};
    if (m == 0) return std::vector<int>(n + 1, 0);
    const int sz = W::size_for(n, "mismatch_counts");  // 読むのは添字 m-1〜n-1 だけなので、巡回の長さは N 以上で足りる
    auto is_w = [&](long long x) { return has_wild && x == wild; };
    std::vector<long long> vals, tv;  // 両方に現れる値（片方にしかない値は一致しない）
    for (long long x : p)
        if (!is_w(x)) vals.push_back(x);
    for (long long x : t)
        if (!is_w(x)) tv.push_back(x);
    std::sort(vals.begin(), vals.end()), std::sort(tv.begin(), tv.end());
    vals.erase(std::unique(vals.begin(), vals.end()), vals.end());
    vals.erase(std::remove_if(vals.begin(), vals.end(), [&](long long x) { return !std::binary_search(tv.begin(), tv.end(), x); }), vals.end());
    std::vector<unsigned> acc(sz, 0), fa(sz), fb(sz);
    auto add_product = [&](auto ft, auto fp, int sign) {  // acc += sign × F(text の印) × F(pat を反転した印)
        std::fill(fa.begin(), fa.end(), 0u), std::fill(fb.begin(), fb.end(), 0u);
        for (int i = 0; i < n; ++i) fa[i] = ft(t[i]) ? 1u : 0u;
        for (int j = 0; j < m; ++j) fb[m - 1 - j] = fp(p[j]) ? 1u : 0u;
        W::ntt(fa), W::ntt(fb);
        for (int k = 0; k < sz; ++k) {
            const unsigned v = (unsigned)((unsigned long long)fa[k] * fb[k] % W::MOD);
            acc[k] = sign > 0 ? (acc[k] + v >= W::MOD ? acc[k] + v - W::MOD : acc[k] + v) : (acc[k] >= v ? acc[k] - v : acc[k] + W::MOD - v);
        }
    };
    // 不一致 = （両方とも wild でない組の数）−（同じ値の組の数）
    if (has_wild) add_product([&](long long x) { return !is_w(x); }, [&](long long x) { return !is_w(x); }, +1);
    for (long long c : vals) add_product([&](long long x) { return x == c; }, [&](long long x) { return x == c; }, -1);
    W::ntt_inv(acc);
    std::vector<int> r(n - m + 1);
    for (int i = 0; i <= n - m; ++i) {
        const long long both = has_wild ? 0 : m;  // wild がなければ組はいつも m 個
        long long v = (long long)acc[i + m - 1];
        if (v > (long long)W::MOD / 2) v -= W::MOD;  // 負の値（−一致の数）を戻す
        r[i] = (int)(both + v);
    }
    return r;
}
template <class T, class P>
std::vector<int> wildcard_match_impl(const T& text0, const P& pat0, long long wild) {
    using W = WcNtt;
    const auto t = codes_of(text0), p = codes_of(pat0);
    const int n = (int)t.size(), m = (int)p.size();
    std::vector<int> res;
    if (m > n) return res;
    if (m == 0) {
        for (int i = 0; i <= n; ++i) res.push_back(i);
        return res;
    }
    const int sz = W::size_for(n, "wildcard_match");
    std::vector<char> ok(n - m + 1, 1);
    static std::mt19937_64 rng((unsigned long long)std::chrono::steady_clock::now().time_since_epoch().count() ^ 0x5DEECE66DULL);
    for (int round = 0; round < 2; ++round) {
        const unsigned long long salt = rng();
        auto rv = [&](long long x) -> unsigned long long {  // 値 → [1, MOD) の乱数（wild は 0）
            if (x == wild) return 0;
            unsigned long long z = (unsigned long long)x ^ salt;
            z += 0x9E3779B97F4A7C15ULL, z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL, z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL, z ^= z >> 31;
            return 1 + z % (W::MOD - 1);
        };
        // Σ_j P_j T_{i+j} (P_j - T_{i+j})^2 = Σ P^3 T − 2 P^2 T^2 + P T^3
        std::vector<unsigned> t1(sz, 0), t2(sz, 0), t3(sz, 0), p1(sz, 0), p2(sz, 0), p3(sz, 0);
        for (int i = 0; i < n; ++i) {
            const unsigned long long v = rv(t[i]);
            t1[i] = (unsigned)v, t2[i] = (unsigned)(v * v % W::MOD), t3[i] = (unsigned)((unsigned long long)t2[i] * v % W::MOD);
        }
        for (int j = 0; j < m; ++j) {
            const unsigned long long v = rv(p[j]);
            const int k = m - 1 - j;
            p1[k] = (unsigned)v, p2[k] = (unsigned)(v * v % W::MOD), p3[k] = (unsigned)((unsigned long long)p2[k] * v % W::MOD);
        }
        W::ntt(t1), W::ntt(t2), W::ntt(t3), W::ntt(p1), W::ntt(p2), W::ntt(p3);
        for (int k = 0; k < sz; ++k) {
            const unsigned long long a = (unsigned long long)p3[k] * t1[k] % W::MOD, b = (unsigned long long)p2[k] * t2[k] % W::MOD,
                                     c = (unsigned long long)p1[k] * t3[k] % W::MOD;
            t1[k] = (unsigned)((a + c + 2ULL * (W::MOD - b)) % W::MOD);
        }
        W::ntt_inv(t1);
        for (int i = 0; i <= n - m; ++i)
            if (t1[i + m - 1] != 0) ok[i] = 0;
    }
    for (int i = 0; i <= n - m; ++i)
        if (ok[i]) res.push_back(i);
    return res;
}
}  // namespace str_internal

// ずらし i ごとの不一致の数（r[i] = pat[j] ≠ text[i + j] の j の数。長さ N − M + 1）  O(σ' N log N)（σ' は text と pat の両方に現れる値の数）
template <class T, class P>
std::vector<int> mismatch_counts(const T& text, const P& pat) { return str_internal::mismatch_counts_impl(text, pat, false, 0); }
// wild つき：両方とも wild でない位置で値が違う数（wild は '?' のような文字か、数の列なら値）
template <class T, class P, class W>
std::vector<int> mismatch_counts(const T& text, const P& pat, W wild) { return str_internal::mismatch_counts_impl(text, pat, true, str_internal::code_of(wild)); }
template <class T, class P>
std::vector<int> match_with_k_mismatches(const T& text, const P& pat, int k) {  // 不一致が k 以下のずらし i の一覧（wild なし）
    std::vector<int> r;
    const auto mc = mismatch_counts(text, pat);
    for (int i = 0; i < (int)mc.size(); ++i)
        if (mc[i] <= k) r.push_back(i);
    return r;
}
template <class T, class P, class W>
std::vector<int> match_with_k_mismatches(const T& text, const P& pat, int k, W wild) {  // wild つき
    std::vector<int> r;
    const auto mc = mismatch_counts(text, pat, wild);
    for (int i = 0; i < (int)mc.size(); ++i)
        if (mc[i] <= k) r.push_back(i);
    return r;
}
// wild を任意の 1 文字として一致するずらし i の一覧（乱数を使う。違うのに一致とする確率は 1 か所あたり約 10^-17）  O(N log N)
template <class T, class P, class W>
std::vector<int> wildcard_match(const T& text, const P& pat, W wild) { return str_internal::wildcard_match_impl(text, pat, str_internal::code_of(wild)); }
template <class T, class P>
std::vector<int> wildcard_match(const T& text, const P& pat) {  // wild = '?'（文字列だけ。数の列では wild の値を渡す）
    static_assert(std::is_same<str_internal::elem_t<T>, char>::value && std::is_same<str_internal::elem_t<P>, char>::value,
                  "数の列では wild の値を渡す：wildcard_match(text, pat, wild)");
    return str_internal::wildcard_match_impl(text, pat, '?');
}
#endif  // GG_WILDCARD_MATCH
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト
// =============================================================================
namespace wc_test {

#define WC_CHECK(...)                                                                          \
    do {                                                                                        \
        if (!(__VA_ARGS__)) {                                                                   \
            std::cerr << "[wildcard_match] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
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

template <class V>
static std::vector<int> bf_mismatch(const V& t, const V& p, bool has_wild, long long wild = 0) {  // 愚直（has_wild = false なら wild なし）
    const int n = (int)t.size(), m = (int)p.size();
    std::vector<int> r;
    if (m > n) return r;
    for (int i = 0; i + m <= n; ++i) {
        int c = 0;
        for (int j = 0; j < m; ++j) {
            const long long a = str_internal::code_of(t[i + j]), b = str_internal::code_of(p[j]);
            if (!(has_wild && (a == wild || b == wild)) && a != b) ++c;
        }
        r.push_back(c);
    }
    return r;
}

static void test_string(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 2000 * TEST_SCALE; ++it) {
        const int sg = 1 + (int)(rng() % 3);
        std::string t = rnd_str_shaped(rng, (int)(rng() % 40), sg), p = rnd_str_shaped(rng, (int)(rng() % 8), sg);
        const int wt = (int)(rng() % 3);  // wild を入れる場所（0 なし・1 pat だけ・2 両方）
        if (wt >= 1)
            for (auto& c : p)
                if (rng() % 4 == 0) c = '?';
        if (wt == 2)
            for (auto& c : t)
                if (rng() % 6 == 0) c = '?';
        // 不一致の数
        const auto bf = bf_mismatch(t, p, true, '?');
        WC_CHECK(mismatch_counts(t, p, '?') == bf);
        const auto bf0 = bf_mismatch(t, p, false);
        WC_CHECK(mismatch_counts(t, p) == bf0);
        const int k = (int)(rng() % 3);
        std::vector<int> ek, ew;
        for (int i = 0; i < (int)bf.size(); ++i) {
            if (bf[i] <= k) ek.push_back(i);
            if (bf[i] == 0) ew.push_back(i);
        }
        WC_CHECK(match_with_k_mismatches(t, p, k, '?') == ek);
        WC_CHECK(wildcard_match(t, p) == ew && wildcard_match(t, p, '?') == ew);
        ++cases;
    }
}

static void test_values(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 1000 * TEST_SCALE; ++it) {
        std::vector<long long> t(rng() % 30), p(rng() % 6);
        const bool big = rng() % 2;
        for (auto& x : t) x = big ? rnd_ll(rng) % 5 * 1000000000000LL : (long long)(rng() % 4) - 1;
        for (auto& x : p) x = big ? rnd_ll(rng) % 5 * 1000000000000LL : (long long)(rng() % 4) - 1;
        if (!big && !p.empty() && rng() % 2) p[rng() % p.size()] = -1;  // -1 を wild に
        WC_CHECK(mismatch_counts(t, p, -1) == bf_mismatch(t, p, true, -1) && mismatch_counts(t, p) == bf_mismatch(t, p, false));
        std::vector<int> ew;
        const auto bf = bf_mismatch(t, p, true, -1);
        for (int i = 0; i < (int)bf.size(); ++i)
            if (bf[i] == 0) ew.push_back(i);
        WC_CHECK(wildcard_match(t, p, -1) == ew);
        ++cases;
    }
    // 端の値（LLONG_MIN・LLONG_MAX）・128 以上の char を wild に・pat にしかない値
    const long long LO = std::numeric_limits<long long>::min(), HI = std::numeric_limits<long long>::max();
    for (int it = 0; it < 500 * TEST_SCALE; ++it) {
        std::vector<long long> t(rng() % 20), p(rng() % 5);
        const long long pool[5] = {LO, HI, 0, -1, (long long)(rng() % 1000) + 5};
        for (auto& x : t) x = pool[rng() % 4];
        for (auto& x : p) x = pool[rng() % 5];
        const long long w = pool[rng() % 4];
        WC_CHECK(mismatch_counts(t, p) == bf_mismatch(t, p, false) && mismatch_counts(t, p, w) == bf_mismatch(t, p, true, w));
        std::vector<int> ew;
        const auto bf = bf_mismatch(t, p, true, w);
        for (int i = 0; i < (int)bf.size(); ++i)
            if (bf[i] == 0) ew.push_back(i);
        WC_CHECK(wildcard_match(t, p, w) == ew && match_with_k_mismatches(t, p, 0, w) == ew);
        std::string a = rnd_str(rng, (int)(rng() % 20), 2), b = rnd_str(rng, (int)(rng() % 5), 3);
        for (auto& c : a)
            if (rng() % 4 == 0) c = '\xff';
        for (auto& c : b)
            if (rng() % 3 == 0) c = '\xff';
        const auto bfs = bf_mismatch(a, b, true, 255);
        std::vector<int> es;
        for (int i = 0; i < (int)bfs.size(); ++i)
            if (bfs[i] == 0) es.push_back(i);
        WC_CHECK(mismatch_counts(a, b, '\xff') == bfs && wildcard_match(a, b, '\xff') == es);
        ++cases;
    }
    // 長い列（NTT の大きさが変わる所）
    for (int it = 0; it < 20 * TEST_SCALE; ++it) {
        const std::string t = rnd_str(rng, 500 + (int)(rng() % 3000), 2), p = rnd_str(rng, 1 + (int)(rng() % 300), 2);
        WC_CHECK(mismatch_counts(t, p) == bf_mismatch(t, p, false));
        ++cases;
    }
}

static void random_test() {
    std::mt19937_64 rng(20261012);
    int cases = 0;
    test_string(rng, cases);
    test_values(rng, cases);
    std::cout << "[wildcard_match] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[wildcard_match] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    auto ms = [](std::chrono::steady_clock::time_point a) {
        return (long long)std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count();
    };
    const std::string s = rnd_str(rng, 1000000, 2, '0'), t = rnd_str(rng, 500000, 2, '0');
    auto t0 = std::chrono::steady_clock::now();
    const auto mc = mismatch_counts(s, t);
    const long long t_mc = ms(t0);
    std::string a = rnd_str(rng, 500000, 26), b = rnd_str(rng, 200000, 26);
    for (auto& c : b)
        if (rng() % 3 == 0) c = '?';
    t0 = std::chrono::steady_clock::now();
    const auto wm = wildcard_match(a, b);
    const long long t_wm = ms(t0);
    t0 = std::chrono::steady_clock::now();
    const auto mc26 = mismatch_counts(a, b, '?');
    const long long t_mc26 = ms(t0);
    std::cout << "[wildcard_match] benchmark: mismatch_counts（σ = 2、10^6 と 5×10^5）" << t_mc << " ms, wildcard_match（σ = 26、5×10^5 と 2×10^5）" << t_wm
              << " ms, mismatch_counts（同じ大きさ・σ = 26）" << t_mc26 << " ms  (check " << (mc[0] + wm.size() + mc26[0]) % 7 << ")\n";
}

}  // namespace wc_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace wc_demo {

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
#define SHOW(expr, comment) wc_demo::show(#expr, (expr), comment)
#define RUN(...)                   \
    wc_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("ワイルドカード・不一致の数（wildcard_match） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  mismatch_counts / match_with_k_mismatches / wildcard_match\n"
              << "  ずらし i は pat を text[i, i + M) に重ねたもの（0〜N - M）\n";

    api("【1】 mismatch_counts(text, pat) / mismatch_counts(text, pat, wild)", "vector<int>（長さ N - M + 1）",
        "ずらしごとに文字が違う数 / wild つき：両方とも wild でない位置で違う数（σ' は text と pat の両方に現れる値の数）", "O(σ' N log N)");
    SHOW(mismatch_counts("0001", "101"), "ABC196 F（Substring 2）の入力例 1 → 最小は 1");
    SHOW(mismatch_counts("abcab", "a?b", '?'), "? はどの文字とも一致");
    {
        const auto mc = mismatch_counts("10101000010011011110", "0010011111");
        SHOW(*std::min_element(mc.begin(), mc.end()), "ABC196 F の入力例 3");
    }

    api("【2】 match_with_k_mismatches(text, pat, k) / (text, pat, k, wild)", "vector<int>", "不一致が k 個以下のずらしの一覧", "同上");
    SHOW(match_with_k_mismatches("abcabd", "abd", 1), "");

    api("【3】 wildcard_match(text, pat) / wildcard_match(text, pat, wild)", "vector<int>",
        "wild（文字列で省くと '?'）を任意の 1 文字として一致するずらしの一覧（両方の wild が使える。乱数を使う。違うのに一致とする確率は 1 か所あたり約 10^-17）", "O(N log N)");
    SHOW(wildcard_match("abcab?b", "a?b"), "text の ? も任意の 1 文字");
    SHOW(wildcard_match(std::vector<int>{3, 1, 4, 1, 5}, std::vector<int>{0, 1}, 0), "数の列では wild の値（ここでは 0）を必ず渡す（省くとコンパイルエラー）");
}

}  // namespace wc_demo

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
    wc_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    wc_test::benchmark();    // 速度計測
    wc_demo::demo();         // 使い方デモ
    return 0;
}
