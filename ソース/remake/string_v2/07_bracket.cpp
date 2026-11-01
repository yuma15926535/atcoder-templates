// =============================================================================
//  括弧列（正しいか・対応・深さ・入れ子の木・対応する組を消す・足りない数・並べ替えてつなげる・区間の判定・書き換え・全列挙）   --- 文字列第二世代
// -----------------------------------------------------------------------------
//  分類       : 07 括弧列
//  依存       : 00_str_base（同梱済み）
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）
//  文字コード : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  括弧の種類は文字列 pairs で渡す（既定 "()"。"()[]<>" なら 3 種類。偶数番目が開き、次がそれに対応する閉じ）
//  括弧以外の文字は無視する（中身として扱う）。対応は「隣り合う同じ種類の開き・閉じを消すことを繰り返す」で一緒に消える組
//  （1 種類）と書いたものは pairs の長さが 2 のときだけ使える（違えば止める）
//
//  機能一覧                                                                       計算量
//    is_balanced(s, pairs)                 正しい括弧列か                                               O(N)
//    match_brackets(s, pairs)              対応する括弧の位置（なければ -1）                              O(N)
//    bracket_depths(s, pairs)              各位置を囲む組の数                                            O(N)
//    bracket_parents(s, pairs)             各位置を囲む一番内側の組の開き括弧（入れ子を木に）               O(N)
//    erase_matched_pairs(s, pairs)         対応する組を中身ごと消した残り（ABC307 D）                     O(N)
//    BracketBalance{need_open, need_close} 対応しない閉じ / 開きの数。+ で連結、balanced()（1 種類）        O(1)
//    balance_of(s, pairs)                  その計算（1 文字 balance_of(c) も）（1 種類）                   O(N)
//    bracket_concat_order(vs, pairs)       並べ替えてつなげて正しくする順 BracketOrder{ok, order}（1 種類）   O(Σ|vs| + K log K)
//    min_insertions_to_balance(s, pairs)   足す括弧の最小の数（1 種類）                                   O(N)
//    fix_brackets(s, pairs)                最小の数を足して正しくしたもの（1 種類）                        O(N)
//    count_balanced_substrings(s, pairs)   正しい括弧列の部分文字列の数                                   O(N)
//    longest_balanced_substring(s, pairs)  最長のもの（StrSpan）                                         O(N)
//    StaticBracket sb(s, pairs)            sb.is_balanced(l, r)（O(1)）/ sb.partner(i)                    前処理 O(N)
//    BracketSegTree bt(s, pairs)           set(i, c) / flip(i) / get(i) / is_balanced(l, r) / balance(l, r) / all_balanced()（1 種類）  O(log N)
//    for_each_balanced(n, f, pairs)        長さ n の正しい括弧列を辞書順（開き < 閉じ）に f(s)（1 種類）    O(個数 × n)
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    brk_test::random_test()  … 隣り合う組を消す愚直解・全部の部分文字列・全部の並べ方と全機能を突き合わせるランダムテスト
//    brk_test::benchmark()    … 速度計測
//    brk_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
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
#ifndef GG_BRACKET
#define GG_BRACKET
// ---- 括弧列（文字列第二世代）  判定・対応・深さ・足りない括弧・つなげる順・区間の判定（書き換えも）・? を決める方法の数・全列挙 ----
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

template <class T, class Seq>
T count_balanced_completions(const Seq& s0, char wild = '?', const std::string& pairs = "()") {  // wild をそれぞれ開き・閉じに決めて正しい括弧列にする方法の数（1 種類。ABC312 D）  O(N^2)
    const auto& s = str_internal::seq(s0);
    const char* func = "count_balanced_completions";
    str_internal::need_single(pairs, func);
    const long long op = (unsigned char)pairs[0], cl = (unsigned char)pairs[1], w = (unsigned char)wild;
    if (w == op || w == cl) str_internal::fail(func, "wild が括弧と同じ文字");
    const int n = (int)s.size();
    std::vector<T> dp(n + 2, T(0)), nd(n + 2, T(0));  // dp[j] = 開き − 閉じ が j の数
    std::vector<char> ov(n + 2, 0), nov(n + 2, 0);
    dp[0] = T(1);
    for (int i = 0; i < n; ++i) {
        const long long c = str_internal::code_of(s[i]);
        if (c != op && c != cl && c != w) str_internal::fail(func, "開き・閉じ・wild 以外の文字がある");
        std::fill(nd.begin(), nd.end(), T(0)), std::fill(nov.begin(), nov.end(), 0);
        for (int j = 0; j <= i; ++j) {
            if (c != cl) {  // 開き
                bool o = nov[j + 1];
                str_internal::add_count(nd[j + 1], o, dp[j], ov[j]);
                nov[j + 1] = o;
            }
            if (c != op && j > 0) {  // 閉じ
                bool o = nov[j - 1];
                str_internal::add_count(nd[j - 1], o, dp[j], ov[j]);
                nov[j - 1] = o;
            }
        }
        std::swap(dp, nd), std::swap(ov, nov);
    }
    if (ov[0]) str_internal::count_overflow(func);
    return dp[0];
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
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト
// =============================================================================
namespace brk_test {

#define BRK_CHECK(...)                                                                          \
    do {                                                                                        \
        if (!(__VA_ARGS__)) {                                                                   \
            std::cerr << "[bracket] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
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

static int kind_of(char c, const std::string& pairs) {  // 0 = 括弧以外、k + 1 = 開き、-(k + 1) = 閉じ
    const auto p = pairs.find(c);
    if (p == std::string::npos) return 0;
    return p % 2 == 0 ? (int)(p / 2) + 1 : -(int)(p / 2) - 1;
}
// 愚直：括弧だけを並べ、隣り合う同じ種類の開き・閉じを消すことを繰り返す
static std::vector<int> bf_match(const std::string& s, const std::string& pairs) {
    std::vector<int> m(s.size(), -1), idx;
    for (int i = 0; i < (int)s.size(); ++i)
        if (kind_of(s[i], pairs) != 0) idx.push_back(i);
    for (bool changed = true; changed;) {
        changed = false;
        for (std::size_t k = 0; k + 1 < idx.size(); ++k) {
            const int a = kind_of(s[idx[k]], pairs), b = kind_of(s[idx[k + 1]], pairs);
            if (a > 0 && b == -a) {
                m[idx[k]] = idx[k + 1], m[idx[k + 1]] = idx[k];
                idx.erase(idx.begin() + k, idx.begin() + k + 2);
                changed = true;
                break;
            }
        }
    }
    return m;
}
static bool bf_balanced(const std::string& s, const std::string& pairs) {
    const auto m = bf_match(s, pairs);
    for (int i = 0; i < (int)s.size(); ++i)
        if (kind_of(s[i], pairs) != 0 && m[i] < 0) return false;
    return true;
}
// 足す括弧の最小の数（区間 DP。何種類でも）
static int bf_min_insert(const std::string& s, const std::string& pairs) {
    const int n = (int)s.size();
    std::vector<std::vector<int>> f(n + 1, std::vector<int>(n + 1, 0));  // f[l][r] = s[l, r)
    for (int len = 1; len <= n; ++len)
        for (int l = 0; l + len <= n; ++l) {
            const int r = l + len, k = kind_of(s[l], pairs);
            if (k == 0) {
                f[l][r] = f[l + 1][r];
                continue;
            }
            int best = 1 + f[l + 1][r];
            if (k > 0)
                for (int m = l + 1; m < r; ++m)
                    if (kind_of(s[m], pairs) == -k) best = std::min(best, f[l + 1][m] + f[m + 1][r]);
            f[l][r] = best;
        }
    return f[0][n];
}
static std::string rnd_brackets(std::mt19937_64& rng, int n, const std::string& pairs, bool letters) {
    std::string s;
    const int t = (int)(rng() % 3);
    if (t == 0) {  // 一様
        for (int i = 0; i < n; ++i) s += (letters && rng() % 4 == 0) ? 'x' : pairs[rng() % pairs.size()];
    } else {  // 正しい括弧列を作ってから少し壊す
        std::string st;
        for (int i = 0; i < n; ++i) {
            if (letters && rng() % 5 == 0) {
                s += 'x';
            } else if (!st.empty() && (rng() % 2 || (int)st.size() >= n - i)) {
                s += pairs[pairs.find(st.back()) + 1], st.pop_back();
            } else {
                const char o = pairs[2 * (rng() % (pairs.size() / 2))];
                s += o, st += o;
            }
        }
        for (int k = (t == 1 ? 0 : 1 + (int)(rng() % 2)); k > 0 && n > 0; --k) s[rng() % n] = pairs[rng() % pairs.size()];
    }
    return s;
}

static void test_basic(std::mt19937_64& rng, int& cases) {
    const std::string kinds[] = {"()", "[]", "()[]", "()[]<>", ")("};
    for (int it = 0; it < 4000 * TEST_SCALE; ++it) {
        const std::string pairs = kinds[rng() % 5];
        const bool letters = rng() % 2;
        const int n = (int)(rng() % (it % 10 == 0 ? 40 : 14));
        const std::string s = rnd_brackets(rng, n, pairs, letters);
        const auto m = bf_match(s, pairs);
        BRK_CHECK(match_brackets(s, pairs) == m);
        BRK_CHECK(is_balanced(s, pairs) == bf_balanced(s, pairs));
        // 深さ・親
        const auto dep = bracket_depths(s, pairs), par = bracket_parents(s, pairs);
        for (int i = 0; i < n; ++i) {
            int d = 0, p = -1;
            for (int l = 0; l < n; ++l)
                if (m[l] > l && l < i && i < m[l]) ++d, p = std::max(p, l);
            BRK_CHECK(dep[i] == d && par[i] == p);
        }
        // 組を消す：繰り返し「開き + 括弧以外 + 閉じ」を消す（ABC307 D の操作）
        {
            std::string t = s;
            for (bool changed = true; changed;) {
                changed = false;
                for (int l = 0; l < (int)t.size() && !changed; ++l) {
                    const int k = kind_of(t[l], pairs);
                    if (k <= 0) continue;
                    int r = l + 1;
                    while (r < (int)t.size() && kind_of(t[r], pairs) == 0) ++r;
                    if (r < (int)t.size() && kind_of(t[r], pairs) == -k) t.erase(l, r - l + 1), changed = true;
                }
            }
            BRK_CHECK(erase_matched_pairs(s, pairs) == t);
        }
        // 部分文字列
        long long cnt = 0;
        StrSpan best{0, 0};
        const StaticBracket sb(s, pairs);
        for (int l = 0; l <= n; ++l)
            for (int r = l; r <= n; ++r) {
                const bool ok = bf_balanced(s.substr(l, r - l), pairs);
                BRK_CHECK(sb.is_balanced(l, r) == ok);
                if (ok && r > l) {
                    ++cnt;
                    if (r - l > best.len()) best = {l, r};
                }
            }
        for (int i = 0; i < n; ++i) BRK_CHECK(sb.partner(i) == m[i]);
        BRK_CHECK(count_balanced_substrings(s, pairs) == cnt);
        BRK_CHECK(longest_balanced_substring(s, pairs) == best);
        // 1 種類
        if (pairs.size() == 2) {
            const BracketBalance b = balance_of(s, pairs);
            long long cur = 0, mn = 0;
            for (char c : s) {
                const int k = kind_of(c, pairs);
                cur += k > 0 ? 1 : (k < 0 ? -1 : 0);
                mn = std::min(mn, cur);
            }
            BRK_CHECK(b.need_open == -mn && b.need_close == cur - mn && b.balanced() == bf_balanced(s, pairs));
            BRK_CHECK(min_insertions_to_balance(s, pairs) == bf_min_insert(s, pairs));
            const std::string fx = fix_brackets(s, pairs);
            BRK_CHECK(bf_balanced(fx, pairs) && (long long)fx.size() == n + min_insertions_to_balance(s, pairs) && is_subsequence(s, fx));
            // 連結
            const int cut = (int)(rng() % (n + 1));
            BRK_CHECK(balance_of(s.substr(0, cut), pairs) + balance_of(s.substr(cut), pairs) == b);
            BracketBalance acc;
            for (char c : s) acc += balance_of(c, pairs);
            BRK_CHECK(acc == b);
        } else {
            BRK_CHECK(bf_min_insert(s, pairs) >= 0);
        }
        ++cases;
    }
}

static void test_concat(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 2000 * TEST_SCALE; ++it) {
        const int k = 1 + (int)(rng() % 6);
        std::vector<std::string> vs(k);
        for (auto& v : vs) v = rnd_brackets(rng, 1 + (int)(rng() % 6), "()", rng() % 4 == 0);
        std::vector<int> perm(k);
        std::iota(perm.begin(), perm.end(), 0);
        bool exists = false;
        do {
            std::string t;
            for (int i : perm) t += vs[i];
            if (bf_balanced(t, "()")) exists = true;
        } while (!exists && std::next_permutation(perm.begin(), perm.end()));
        const BracketOrder r = bracket_concat_order(vs);
        BRK_CHECK(r.ok == exists);
        if (r.ok) {
            BRK_CHECK((int)r.order.size() == k);
            std::vector<int> seen(k, 0);
            std::string t;
            for (int i : r.order) {
                BRK_CHECK(i >= 0 && i < k && !seen[i]);
                seen[i] = 1, t += vs[i];
            }
            BRK_CHECK(bf_balanced(t, "()"));
        } else {
            BRK_CHECK(r.order.empty());
        }
        ++cases;
    }
}

static void test_segtree(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 300 * TEST_SCALE; ++it) {
        const std::string pairs = rng() % 2 ? "()" : "<>";
        const int n = 1 + (int)(rng() % 40);
        std::string s = rnd_brackets(rng, n, pairs, rng() % 2);
        BracketSegTree bt(s, pairs);
        for (int q = 0; q < 200; ++q) {
            const int op = (int)(rng() % 4);
            if (op == 0) {
                const int i = (int)(rng() % n);
                const char c = rng() % 5 == 0 ? 'x' : pairs[rng() % 2];
                s[i] = c, bt.set(i, c);
            } else if (op == 1) {
                const int i = (int)(rng() % n);
                if (s[i] == pairs[0]) s[i] = pairs[1];
                else if (s[i] == pairs[1]) s[i] = pairs[0];
                bt.flip(i);
            } else {
                int l = (int)(rng() % (n + 1)), r = (int)(rng() % (n + 1));
                if (l > r) std::swap(l, r);
                const std::string sub = s.substr(l, r - l);
                BRK_CHECK(bt.is_balanced(l, r) == bf_balanced(sub, pairs) && bt.balance(l, r) == balance_of(sub, pairs));
            }
            const int i = (int)(rng() % n);
            BRK_CHECK(bt.get(i) == s[i] && bt.all_balanced() == bf_balanced(s, pairs) && bt.to_string() == s);
        }
        ++cases;
    }
}

static void test_enumerate(int& cases) {
    for (const std::string pairs : {"()", ")(", "[]"}) {
        for (int n = 0; n <= 14; ++n) {
            std::vector<std::string> got, want;
            for_each_balanced(n, [&](const std::string& s) { got.push_back(s); }, pairs);
            for (int mask = 0; mask < (1 << n); ++mask) {
                std::string s(n, pairs[0]);
                for (int i = 0; i < n; ++i)
                    if (mask >> (n - 1 - i) & 1) s[i] = pairs[1];  // mask の昇順 = 開き < 閉じ の辞書順
                if (bf_balanced(s, pairs)) want.push_back(s);
            }
            BRK_CHECK(got == want);
            ++cases;
        }
    }
}

static void test_final_additions(std::mt19937_64& rng, int& cases) {  // count_balanced_completions
    for (int it = 0; it < 2000 * TEST_SCALE; ++it) {
        const int n = (int)(rng() % 13);
        std::string s(n, ' ');
        for (char& c : s) c = "()?"[rng() % 3];
        std::vector<int> q;
        for (int i = 0; i < n; ++i)
            if (s[i] == '?') q.push_back(i);
        long long cnt = 0;
        for (int mask = 0; mask < (1 << q.size()); ++mask) {
            std::string t = s;
            for (std::size_t k = 0; k < q.size(); ++k) t[q[k]] = (mask >> k & 1) ? '(' : ')';
            cnt += is_balanced(t);
        }
        BRK_CHECK(count_balanced_completions<long long>(s) == cnt && count_balanced_completions<double>(s) == (double)cnt);
        std::string u = s;
        for (char& c : u) c = c == '(' ? 'A' : c == ')' ? 'B' : '*';
        BRK_CHECK(count_balanced_completions<long long>(u, '*', "AB") == cnt);
        ++cases;
    }
    // 途中だけ溢れても答えが収まれば止まらない：? が 70 個で答えは C(70, 35)/36（カタラン数）
    BRK_CHECK(count_balanced_completions<long long>(std::string(70, '?')) == 3116285494907301262LL);
}
static void random_test() {
    std::mt19937_64 rng(20261012);
    int cases = 0;
    test_basic(rng, cases);
    test_concat(rng, cases);
    test_segtree(rng, cases);
    test_final_additions(rng, cases);
    test_enumerate(cases);
    std::cout << "[bracket] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[bracket] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    auto ms = [](std::chrono::steady_clock::time_point a) {
        return (long long)std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count();
    };
    const int N = 1000000;
    const std::string s = rnd_brackets(rng, N, "()[]", false);
    auto t0 = std::chrono::steady_clock::now();
    long long x = is_balanced(s, "()[]");
    x += match_brackets(s, "()[]")[N / 2];
    x += count_balanced_substrings(s, "()[]");
    const long long t_basic = ms(t0);
    t0 = std::chrono::steady_clock::now();
    const StaticBracket sb(s, "()[]");
    for (int q = 0; q < N; ++q) {
        int l = (int)(rng() % (N + 1)), r = (int)(rng() % (N + 1));
        if (l > r) std::swap(l, r);
        x += sb.is_balanced(l, r);
    }
    const long long t_static = ms(t0);
    const int M = 500000;
    BracketSegTree bt(rnd_brackets(rng, M, "()", false));
    t0 = std::chrono::steady_clock::now();
    for (int q = 0; q < M; ++q) {
        if (q % 2) bt.flip((int)(rng() % M));
        else {
            int l = (int)(rng() % (M + 1)), r = (int)(rng() % (M + 1));
            if (l > r) std::swap(l, r);
            x += bt.is_balanced(l, r);
        }
    }
    const long long t_seg = ms(t0);
    t0 = std::chrono::steady_clock::now();
    long long cnt = 0;
    for_each_balanced(24, [&](const std::string& t) { cnt += t[1] == '('; });
    const long long t_enum = ms(t0);
    std::cout << "[bracket] benchmark: N = 10^6 で is_balanced + match_brackets + count_balanced_substrings " << t_basic << " ms, StaticBracket 作成 + 区間 10^6 回 "
              << t_static << " ms, BracketSegTree（N = Q = 5×10^5）" << t_seg << " ms, for_each_balanced(24)（208012 個）" << t_enum << " ms  (check " << (x + cnt) % 7 << ")\n";
}

}  // namespace brk_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace brk_demo {

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
#define SHOW(expr, comment) brk_demo::show(#expr, (expr), comment)
#define RUN(...)                   \
    brk_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("括弧列（bracket） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  is_balanced / match_brackets / bracket_depths / bracket_parents / erase_matched_pairs   判定・対応・深さ・入れ子の木・組を消す\n"
              << "  BracketBalance / balance_of / bracket_concat_order / min_insertions_to_balance / fix_brackets   足りない括弧の数（1 種類）\n"
              << "  count_balanced_substrings / longest_balanced_substring   正しい部分文字列の数・最長\n"
              << "  StaticBracket（区間の判定 O(1)）/ BracketSegTree（書き換えと区間の判定。1 種類）/ for_each_balanced（全列挙）\n"
              << "  count_balanced_completions<T>（? を決めて正しい括弧列にする方法の数）\n"
              << "  括弧の種類は pairs（既定 \"()\"。\"()[]<>\" なら 3 種類）。括弧以外の文字は無視する\n";

    api("【1】 is_balanced(s, pairs)", "bool", "正しい括弧列か（隣り合う同じ種類の開き・閉じを消して空になるか）", "O(N)");
    SHOW(is_balanced("(()())"), "");
    SHOW(is_balanced("())("), "");
    SHOW(is_balanced("a(b)c"), "括弧以外は無視");
    SHOW(is_balanced("([])<>()", "()[]<>"), "ABC394 D（Colorful Bracket Sequence）の入力例 1");
    SHOW(is_balanced("([<)]>", "()[]<>"), "入力例 2：種類が交差");

    api("【2】 match_brackets(s, pairs) / bracket_depths(s, pairs) / bracket_parents(s, pairs)", "vector<int>",
        "対応する括弧の位置（なければ -1）/ 各位置を囲む組の数 / 囲む一番内側の組の開き括弧の位置（なければ -1）", "O(N)");
    SHOW(match_brackets("(()))("), "最後の 2 つは対応しない");
    SHOW(bracket_depths("(a(b))c"), "");
    SHOW(bracket_parents("(a(b))c"), "");

    api("【3】 erase_matched_pairs(s, pairs)", "string", "対応する括弧の組を中身ごと消した残り", "O(N)");
    SHOW(erase_matched_pairs("a(b(d))c"), "ABC307 D（Mismatched Parentheses）の入力例 1");
    SHOW(erase_matched_pairs("a(b)("), "入力例 2");

    api("【4】 BracketBalance{need_open, need_close} / balance_of(s, pairs)（1 種類）", "BracketBalance",
        "対応しない閉じ括弧の数（前に足す開き）と開き括弧の数（後ろに足す閉じ）。a + b でつなげたもの。balanced() で正しいか", "O(N) / 連結 O(1)");
    SHOW(balance_of("))(()("), "{前に 2 個, 後ろに 2 個}");
    SHOW(balance_of(")(") + balance_of("()"), "");
    SHOW((balance_of("((") + balance_of(")") + balance_of(")")).balanced(), "");
    {
        // ABC428 C：1 文字足す・最後の 1 文字を消す を繰り返して、毎回 正しい括弧列か
        note("ABC428 C（Brackets Stack Query）: 「先頭からの BracketBalance」を積み、足すときは 1 文字ぶんを + して push、消すときは pop");
        std::vector<BracketBalance> st = {BracketBalance{}};
        std::string out;
        const std::vector<std::string> qs = {"1 (", "2", "1 (", "1 )", "2", "1 (", "1 )", "1 )"};
        for (const auto& q : qs) {
            if (q[0] == '1') st.push_back(st.back() + balance_of(q[2]));
            else st.pop_back();
            out += st.back().balanced() ? "Yes " : "No ";
        }
        note("入力例 1 → " + out);
    }

    api("【5】 bracket_concat_order(vs, pairs)（1 種類）", "BracketOrder{ok, order}",
        "全部の文字列を並べ替えてつなげ、正しい括弧列にする順（開きが多いものを必要な開きの少ない順 → 閉じが多いものを残す開きの多い順）", "O(Σ|vs| + K log K)");
    {
        const std::vector<std::string> vs = {")", "(()"};
        const BracketOrder r = bracket_concat_order(vs);
        code("bracket_concat_order({\")\", \"(()\"})   // ABC167 F（Bracket Sequencing）の入力例 1");
        std::cout << "    -> ok = " << to_s(r.ok) << ", order = " << to_s(r.order) << "\n";
        SHOW(bracket_concat_order(std::vector<std::string>{")(", "()"}).ok, "入力例 2");
    }

    api("【6】 min_insertions_to_balance(s, pairs) / fix_brackets(s, pairs)（1 種類）", "long long / string",
        "正しくするのに足す括弧の最小の数 / 前に開き・後ろに閉じを足して正しくしたもの", "O(N)");
    SHOW(min_insertions_to_balance("())("), "");
    SHOW(fix_brackets("())("), "");

    api("【7】 count_balanced_substrings(s, pairs) / longest_balanced_substring(s, pairs)", "long long / StrSpan",
        "正しい括弧列である空でない部分文字列の数 / 最長のもの（同じ長さなら左）", "O(N)");
    SHOW(count_balanced_substrings("()()(()"), "() が 3 つ・()() が 1 つ");
    SHOW(longest_balanced_substring(")()())"), "");

    api("【8】 StaticBracket sb(s, pairs)", "sb.is_balanced(l, r) -> bool / sb.partner(i) -> int",
        "書き換えのない文字列で、区間 [l, r) が正しい括弧列か O(1)（何種類でも）。L 文字目〜R 文字目は (L - 1, R)", "前処理 O(N)・質問 O(1)");
    {
        RUN(const StaticBracket sb("(())("));
        SHOW(sb.is_balanced(0, 4), "ABC223 F の入力例 1：1〜4 文字目");
        SHOW(sb.is_balanced(0, 2), "1〜2 文字目");
        SHOW(sb.partner(1), "");
    }

    api("【9】 BracketSegTree bt(s, pairs)（1 種類）", "set(i, c) / flip(i) / get(i) / is_balanced(l, r) / balance(l, r) / all_balanced()",
        "1 文字の書き換え・開き閉じの反転と、区間が正しい括弧列か（セグメント木）", "O(log N)（all_balanced は O(1)）");
    {
        RUN(BracketSegTree bt("(())("));
        SHOW(bt.is_balanced(0, 4), "ABC223 F の入力例 2：判定");
        note("1 4 の入れ替え: a = bt.get(0), b = bt.get(3) を読んで bt.set(0, b), bt.set(3, a)");
        const char a = bt.get(0), b = bt.get(3);
        bt.set(0, b), bt.set(3, a);
        SHOW(bt.to_string(), "");
        SHOW(bt.is_balanced(0, 4), "");
        RUN(bt.flip(4));
        SHOW(bt.balance(0, 5), "");
    }

    api("【10】 for_each_balanced(n, f, pairs)（1 種類）", "なし（f(const string&) を呼ぶ）", "長さ n の正しい括弧列を辞書順（開き < 閉じ）に全部。n が奇数なら呼ばない", "O(個数 × n)");
    {
        std::vector<std::string> all;
        for_each_balanced(6, [&](const std::string& t) { all.push_back(t); });
        note("典型90 002（Encyclopedia of Parentheses）の形。n = 6 → " + join(all));
    }

    head("問題での使い方：ABC350 F（Transpose）");
    note("組の中身を反転して大文字小文字を入れ替える → m = match_brackets(S)、d = bracket_depths(S) を作り、");
    note("i = 0 から向き +1 で歩く。括弧なら i = m[i] に飛んで向きを反転、文字なら d[i] が奇数なら入れ替えて出力。");
    {
        const std::string S = "((XYZ)n(X(y)Z))";
        const auto m = match_brackets(S);
        const auto d = bracket_depths(S);
        std::string out;
        for (int i = 0, dir = 1; 0 <= i && i < (int)S.size(); i += dir) {
            if (m[i] >= 0) i = m[i], dir = -dir;
            else out += d[i] % 2 ? swap_case(std::string(1, S[i]))[0] : S[i];
        }
        note("入力例 2 \"" + S + "\" → " + out);
    }

    api("【11】 count_balanced_completions<T>(s, wild = '?', pairs)（1 種類）", "T（mint・long long。答えが溢れたら止める）",
        "wild をそれぞれ開き・閉じに決めて正しい括弧列にする方法の数（dp[開き − 閉じ]）", "O(N^2)");
    SHOW(count_balanced_completions<long long>(std::string("(??") + "?(?"), "ABC312 D の入力例 1（2 つに分けて書いたのは ? の 3 連続の警告を避けるため）");
    SHOW(count_balanced_completions<long long>(std::string(10, '?')), "カタラン数 C_5");
    note("答えを mod で出すときは count_balanced_completions<mint>(S)。");

    head("【型】 途中で閉じが開きより多くならないか（区間・書き換えあり。ABC473 F）");
    note("「正しい括弧列」ではなく「どの接頭辞でも 閉じ ≤ 開き（最後に開きが余ってよい）」は BracketSegTree の balance(l, r).need_open == 0。");
    note("括弧でない 2 文字も pairs に渡せる：BracketSegTree bt(S, \"AB\")（A が開き、B が閉じ）。");
    {
        BracketSegTree bt("AABAB", "AB");
        std::cout << "    S = AABAB：balance(0, 5).need_open = " << bt.balance(0, 5).need_open << "（途中で B が多くならない）、balance(2, 5).need_open = "
                  << bt.balance(2, 5).need_open << "（\"BAB\" は先頭の B が多い）\n";
    }
}

}  // namespace brk_demo

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
    brk_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    brk_test::benchmark();    // 速度計測
    brk_demo::demo();         // 使い方デモ
    return 0;
}
