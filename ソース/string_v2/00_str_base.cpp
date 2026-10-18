// =============================================================================
//  文字列基盤（文字の種類・区間の文字の数・出現位置・ランレングス・変換・比較・区切り）   --- 文字列第二世代
// -----------------------------------------------------------------------------
//  分類       : 00 基盤（ほかの文字列テンプレートはすべてこれを同梱している）
//  依存       : なし（ほかの第二世代（グラフ・グリッド・幾何・数学）と一緒に貼っても衝突しない）
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）
//  文字コード : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  共通の約束（第二世代の文字列テンプレートすべてで同じ）
//    - 列を受け取る関数・構造体は string と vector<整数>（int・long long）のどちらでもよい。文字列リテラル "abc" も渡せる
//    - 範囲は半開区間 [l, r)（0 始まり）。問題の「L 文字目から R 文字目」は (L - 1, R) で渡す。k 番目は 0 から
//    - 位置・番号が「ない」ときは -1
//    - 文字の種類は Alphabet で指定する（既定は英小文字 a〜z）。番号の順が辞書順として扱われる
//    - 範囲外の文字・使えない引数（範囲外の位置など）はメッセージを出して止める（黙って配列の外を読まない）
//    - char は unsigned char として扱う（128 以上の文字（日本語など）でも負の添字にしない。std::string の比較と同じ順）
//    - 「部分文字列」は連続する部分、「部分列」は飛び飛びに選んだもの（名前で区別する）
//
//  機能一覧                                                                       計算量
//   文字の種類
//    Alphabet::lower() / upper() / digits()       a〜z（既定）/ A〜Z / 0〜9
//    Alphabet::letters() / alnum()                A〜Z,a〜z（52 種）/ 0〜9,A〜Z,a〜z（62 種）。どちらも ASCII 順
//    Alphabet::ascii() / bytes()                  文字コード 0〜127 / 0〜255
//    Alphabet::of("()[]")                         並べた文字だけ（並べた順が番号 = 辞書順の順。独自の順にも）
//    Alphabet::ints(K) / ints(lo, hi)             整数 0〜K-1 / [lo, hi)（数の列に）
//    A.id(c) / A.contains(c) / A.size()           文字・値の番号（範囲外は止める）/ 入っているか / 種類の数     O(1)
//    A.value(i) / A.ch(i) / A.chars()             番号 → 値 / 文字 / 全部の文字
//    to_ids(s, A) / from_ids(ids, A)              文字列 → 番号の列 / 番号の列 → 文字列                     O(N)
//    compress_values(a)                           数の列を 0〜K-1 に圧縮（CompressedValues{ids, values}。.alphabet() で ints(K)）  O(N log N)
//   文字の数・位置
//    count_chars(s, A)                            番号ごとの個数                                            O(N)
//    CharPrefix cp(s, A)                          区間の文字の数・次 / 前の出現を O(1) で（前処理 O(Nσ)・メモリ (N+1)σ 個の int）
//      cp.count(l, r, c) / has(l, r, c)           [l, r) にある c の数 / あるか                              O(1)
//      cp.counts(l, r) / distinct(l, r)           [l, r) の番号ごとの数 / 種類の数                           O(σ)
//      cp.next(i, c) / prev(i, c)                 i 以上で最初の c / i より前で最後の c（なければ -1）         O(1)
//      cp.kth(c, k) / total(c)                    c の k 番目（0 から）の位置 / 個数                          O(1)
//    CharPositions cps(s, A)                      文字の出現位置（メモリを抑えた版。前処理 O(N + σ)・メモリ O(N + σ)。σ が大きくてもよい）
//      cps.next(i, c) / prev(i, c)                i 以上で最初の c / i より前で最後の c（なければ -1）         O(log N)
//      cps.kth(c, k) / total(c) / positions(c)    c の k 番目（0 から）の位置 / 個数 / 位置の一覧              O(1)
//      cps.count(l, r, c) / has(l, r, c)          [l, r) にある c の数 / あるか                              O(log N)
//      cps.occurrence_index(i)                    s[i] がその文字の何番目（0 から）の出現か                    O(1)
//   ランレングス
//    run_length(s)                                連続する同じ値を RunBlock{value, len, start} の列に        O(N)
//    run_length_decode(runs)                      元に戻す（char なら string、それ以外は vector）
//   変換（元の列は変えず、新しい列を返す）
//    to_lower(s) / to_upper(s) / swap_case(s)     英字の大文字・小文字                                       O(N)
//    caesar_shift(s, k, A)                        A の文字を番号で k ずらす（負も。A にない文字はそのまま）      O(N)
//    CharMap cm(A)                                文字の置き換えを重ねる
//      cm.replace(c, d)                           今 c になっている文字を全部 d に（ABC342 C）                 O(σ)
//      cm.set(c, d) / cm(c) / cm.apply(s)         元の c の行き先を d に / 元の c の行き先 / s に当てはめる     O(1) / O(1) / O(N)
//      cm.compose(g) / cm.table()                 「cm のあとに g」/ 番号 → 行き先の番号
//    reversed(s) / min_with_reverse(s)            反転 / s と反転の小さい方（向きを区別しないときの代表）       O(N)
//    sorted_chars(s) / repeat_string(s, k)        並べ替えたもの / k 回つなげたもの
//   比較
//    common_prefix_length(a, b) / common_suffix_length(a, b)   共通の接頭辞 / 接尾辞の長さ                   O(N)
//    hamming_distance(a, b)                       同じ位置で違う文字の数（同じ長さ）                          O(N)
//    is_subsequence(t, s)                         t が s の部分列か（何度も聞くなら 06_subsequence）           O(|s|)
//    is_palindrome(s)                             全体が回文か（区間を何度も聞くなら 04_palindrome）           O(N)
//   区切り・位置
//    split(s, sep, skip_empty) / join(vs, sep)    区切る / 区切りを挟んでつなげる（数も。答えの出力に）         O(N)
//    StrSpan{l, r} / substr_of(s, sp)             部分文字列の位置（len()）/ その部分のコピー
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    strb_test::random_test()  … 愚直解と全機能を突き合わせるランダムテスト
//    strb_test::benchmark()    … 速度計測
//    strb_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
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
#include <type_traits>
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
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト
// =============================================================================
namespace strb_test {

#define STRB_CHECK(...)                                                                          \
    do {                                                                                        \
        if (!(__VA_ARGS__)) {                                                                   \
            std::cerr << "[str_base] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
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

__extension__ typedef __int128 I128;
static const long long LMAX = std::numeric_limits<long long>::max();

// 愚直な Alphabet：期待する「値 → 番号」を map で持つ
static void check_alphabet(const Alphabet& A, const std::map<long long, int>& expect, int& cases) {
    STRB_CHECK(A.size() == (int)expect.size());
    for (int c = 0; c < 256; ++c) {
        const auto it = expect.find(c);
        STRB_CHECK(A.contains((char)c) == (it != expect.end()));
        STRB_CHECK(A.contains(c) == (it != expect.end()));
        if (it != expect.end()) STRB_CHECK(A.id((char)c) == it->second && A.id(c) == it->second);
        ++cases;
    }
    for (const auto& [v, i] : expect) {
        STRB_CHECK(A.id(v) == i);
        STRB_CHECK(A.value(i) == v);
        if (0 <= v && v < 256) STRB_CHECK((unsigned char)A.ch(i) == v);
    }
}

static void test_alphabet(std::mt19937_64& rng, int& cases) {
    auto chars_map = [](const std::string& cs) {
        std::map<long long, int> m;
        for (int i = 0; i < (int)cs.size(); ++i) m[(unsigned char)cs[i]] = i;
        return m;
    };
    std::string lo, up, dg;
    for (char c = 'a'; c <= 'z'; ++c) lo += c;
    for (char c = 'A'; c <= 'Z'; ++c) up += c;
    for (char c = '0'; c <= '9'; ++c) dg += c;
    check_alphabet(Alphabet(), chars_map(lo), cases);
    check_alphabet(Alphabet::lower(), chars_map(lo), cases);
    check_alphabet(Alphabet::upper(), chars_map(up), cases);
    check_alphabet(Alphabet::digits(), chars_map(dg), cases);
    check_alphabet(Alphabet::letters(), chars_map(up + lo), cases);
    check_alphabet(Alphabet::alnum(), chars_map(dg + up + lo), cases);
    std::string as, bs;
    for (int c = 0; c < 128; ++c) as += (char)c;
    for (int c = 0; c < 256; ++c) bs += (char)c;
    check_alphabet(Alphabet::ascii(), chars_map(as), cases);
    check_alphabet(Alphabet::bytes(), chars_map(bs), cases);
    STRB_CHECK(Alphabet::lower().chars() == lo && Alphabet::bytes().chars() == bs);
    STRB_CHECK(Alphabet::lower().is_chars() && !Alphabet::ints(5).is_chars());
    // letters / alnum / bytes の番号の順は std::string の比較の順と同じ
    for (const Alphabet& A : {Alphabet::letters(), Alphabet::alnum(), Alphabet::bytes()})
        for (int i = 0; i + 1 < A.size(); ++i) STRB_CHECK(std::string(1, A.ch(i)) < std::string(1, A.ch(i + 1)));
    // of：ランダムな文字の並び
    for (int it = 0; it < 300 * TEST_SCALE; ++it) {
        std::string cs = bs;
        std::shuffle(cs.begin(), cs.end(), rng);
        cs.resize(1 + rng() % 256);
        check_alphabet(Alphabet::of(cs), chars_map(cs), cases);
        STRB_CHECK(Alphabet::of(cs).chars() == cs);
    }
    // ints：0〜K-1 と [lo, hi)（大きな値・負の値も）
    for (int it = 0; it < 3000 * TEST_SCALE; ++it) {
        const long long w = 1 + (long long)(rng() % (rng() % 2 ? 50 : 1000000000ULL));
        long long l = rnd_ll(rng);
        if ((I128)l + w > (I128)LMAX) l = LMAX - w;
        const Alphabet A = (it % 3 == 0) ? Alphabet::ints(w) : Alphabet::ints(l, l + w);
        const long long base = (it % 3 == 0) ? 0 : l;
        STRB_CHECK(A.size() == (int)w);
        for (int t = 0; t < 20; ++t) {
            long long v = rng() % 2 ? base + (long long)(rng() % (unsigned long long)w) : rnd_ll(rng);
            if (rng() % 4 == 0) {  // 境目の外側（long long に収まるときだけ）
                const I128 u = (I128)base + (rng() % 2 ? -1 : w);
                if (u >= -(I128)LMAX - 1 && u <= (I128)LMAX) v = (long long)u;
            }
            const bool in = (I128)v >= base && (I128)v < (I128)base + w;
            STRB_CHECK(A.contains(v) == in);
            if (in) STRB_CHECK(A.id(v) == (int)(v - base) && A.value(A.id(v)) == v);
            ++cases;
        }
    }
    // to_ids / from_ids
    for (int it = 0; it < 200 * TEST_SCALE; ++it) {
        const std::string s = rnd_str_shaped(rng, (int)(rng() % 30), 26);
        const auto ids = to_ids(s);
        STRB_CHECK(ids.size() == s.size());
        for (std::size_t i = 0; i < s.size(); ++i) STRB_CHECK(ids[i] == s[i] - 'a');
        STRB_CHECK(from_ids(ids) == s);
        const std::string u = to_upper(s);
        STRB_CHECK(from_ids(to_ids(u, Alphabet::upper()), Alphabet::upper()) == u);
        ++cases;
    }
}

static void test_compress(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 2000 * TEST_SCALE; ++it) {
        const int n = (int)(rng() % 40);
        std::vector<long long> a(n);
        for (auto& x : a) x = rng() % 3 ? (long long)(rng() % 7) - 3 : rnd_ll(rng);
        const CompressedValues c = compress_values(a);
        STRB_CHECK((int)c.ids.size() == n);
        for (int i = 0; i + 1 < c.size(); ++i) STRB_CHECK(c.values[i] < c.values[i + 1]);
        std::set<long long> st(a.begin(), a.end());
        STRB_CHECK((int)st.size() == c.size());
        for (int i = 0; i < n; ++i) STRB_CHECK(c.values[c.ids[i]] == a[i]);
        STRB_CHECK(c.alphabet().size() == std::max(1, c.size()));
        if (n) STRB_CHECK(to_ids(c.ids, c.alphabet()) == c.ids);
        ++cases;
    }
    // string も（文字コードで並ぶ）
    const CompressedValues c = compress_values("banana");
    STRB_CHECK(c.values == std::vector<long long>({'a', 'b', 'n'}) && c.ids == std::vector<int>({1, 0, 2, 0, 2, 0}));
}

// 区間の文字の数・出現位置を愚直と比べる（s は列、A は文字の種類、vals は A の値の一覧）
template <class Seq>
static void check_char_queries(std::mt19937_64& rng, const Seq& s, const Alphabet& A, int& cases) {
    const int n = (int)s.size(), sg = A.size();
    const CharPrefix cp(s, A);
    const CharPositions cps(s, A);
    STRB_CHECK(cp.size() == n && cps.size() == n);
    const auto cnt_all = count_chars(s, A);
    for (int k = 0; k < sg; ++k) {
        int c = 0;
        std::vector<int> pos;
        for (int i = 0; i < n; ++i)
            if (A.id(s[i]) == k) ++c, pos.push_back(i);
        const long long v = A.value(k);
        STRB_CHECK(cnt_all[k] == c && cps.total(v) == c && cp.total(v) == c && cp.total_id(k) == c && cps.positions(v) == pos);
        for (int j = 0; j < (int)pos.size(); ++j) STRB_CHECK(cps.kth(v, j) == pos[j] && cp.kth(v, j) == pos[j] && cp.kth_id(k, j) == pos[j] && cps.occurrence_index(pos[j]) == j);
        STRB_CHECK(cps.kth(v, (long long)pos.size()) == -1 && cps.kth(v, -1) == -1 && cp.kth(v, (long long)pos.size()) == -1 && cp.kth(v, -1) == -1);
    }
    const int Q = n <= 12 ? (n + 1) * (n + 1) : 200;
    for (int q = 0; q < Q; ++q) {
        int l, r;
        if (n <= 12) l = q / (n + 1), r = q % (n + 1);
        else l = (int)(rng() % (n + 1)), r = (int)(rng() % (n + 1));
        if (l > r) continue;
        std::vector<long long> bf(sg, 0);
        for (int i = l; i < r; ++i) ++bf[A.id(s[i])];
        STRB_CHECK(cp.counts(l, r) == bf);
        int d = 0;
        for (int k = 0; k < sg; ++k) d += bf[k] > 0;
        STRB_CHECK(cp.distinct(l, r) == d);
        for (int t = 0; t < 4; ++t) {
            const int k = (int)(rng() % sg);
            const long long v = A.value(k);
            STRB_CHECK(cp.count(l, r, v) == bf[k] && cp.count_id(l, r, k) == bf[k] && cps.count(l, r, v) == bf[k]);
            STRB_CHECK(cp.has(l, r, v) == (bf[k] > 0) && cps.has(l, r, v) == (bf[k] > 0));
        }
        ++cases;
    }
    for (int i = 0; i <= n; ++i) {
        for (int t = 0; t < 3; ++t) {
            const int k = (int)(rng() % sg);
            const long long v = A.value(k);
            int nx = -1, pv = -1;
            for (int j = i; j < n; ++j)
                if (A.id(s[j]) == k) { nx = j; break; }
            for (int j = i - 1; j >= 0; --j)
                if (A.id(s[j]) == k) { pv = j; break; }
            STRB_CHECK(cps.next(i, v) == nx && cps.prev(i, v) == pv && cp.next(i, v) == nx && cp.prev(i, v) == pv && cp.next_id(i, k) == nx && cp.prev_id(i, k) == pv);
        }
    }
}

static void test_char_queries(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 300 * TEST_SCALE; ++it) {
        const int n = (int)(rng() % (it % 4 == 0 ? 60 : 13));
        const int sg = 1 + (int)(rng() % 4);
        // 英小文字（σ を小さく使う）
        check_char_queries(rng, rnd_str_shaped(rng, n, sg), Alphabet::lower(), cases);
        // 括弧（並べた文字）
        std::string br = rnd_str(rng, n, 4);
        for (char& c : br) c = "()[]"[c - 'a'];
        check_char_queries(rng, br, Alphabet::of("()[]"), cases);
        // 128 以上の文字も（bytes）
        std::string by(n, '\0');
        for (char& c : by) c = (char)(unsigned char)(rng() % 2 ? 200 + rng() % 3 : 'a' + rng() % 2);
        check_char_queries(rng, by, Alphabet::bytes(), cases);
        // 整数の列（負の範囲の ints(lo, hi)）
        std::vector<long long> a(n);
        const long long lo = -(long long)(rng() % 5) - 1000000000000LL;
        for (auto& x : a) x = lo + (long long)(rng() % (unsigned)sg);
        check_char_queries(rng, a, Alphabet::ints(lo, lo + sg), cases);
        std::vector<int> b(n);
        for (auto& x : b) x = (int)(rng() % (unsigned)sg);
        check_char_queries(rng, b, Alphabet::ints(sg), cases);
    }
}

static void test_run_length(std::mt19937_64& rng, int& cases) {
    STRB_CHECK(run_length(std::string()).empty() && run_length_decode(run_length(std::string())).empty());
    for (int it = 0; it < 3000 * TEST_SCALE; ++it) {
        const int n = (int)(rng() % 40);
        const std::string s = rnd_str_shaped(rng, n, 1 + (int)(rng() % 3));
        const auto runs = run_length(s);
        STRB_CHECK(run_length_decode(runs) == s);
        long long pos = 0;
        for (std::size_t i = 0; i < runs.size(); ++i) {
            STRB_CHECK(runs[i].start == pos && runs[i].len >= 1);
            if (i) STRB_CHECK(runs[i].value != runs[i - 1].value);
            for (long long j = 0; j < runs[i].len; ++j) STRB_CHECK(s[pos + j] == runs[i].value);
            pos += runs[i].len;
        }
        STRB_CHECK(pos == n);
        std::vector<long long> a(n);
        for (int i = 0; i < n; ++i) a[i] = (s[i] - 'a') * 1000000000000LL - 5;
        const auto ra = run_length(a);
        STRB_CHECK(run_length_decode(ra) == a && ra.size() == runs.size());
        ++cases;
    }
    // 自分で作った長い run を戻す
    std::vector<RunBlock<char>> rb = {{'x', 3, 0}, {'y', 0, 3}, {'z', 2, 3}};
    STRB_CHECK(run_length_decode(rb) == "xxxzz");
}

static void test_transform(std::mt19937_64& rng, int& cases) {
    std::string all;
    for (int c = 0; c < 256; ++c) all += (char)c;
    const std::string lo = to_lower(all), up = to_upper(all), sw = swap_case(all);
    for (int c = 0; c < 256; ++c) {
        const bool U = 'A' <= c && c <= 'Z', L = 'a' <= c && c <= 'z';
        STRB_CHECK((unsigned char)lo[c] == (U ? c + 32 : c));
        STRB_CHECK((unsigned char)up[c] == (L ? c - 32 : c));
        STRB_CHECK((unsigned char)sw[c] == (U ? c + 32 : L ? c - 32 : c));
        STRB_CHECK(to_lower((char)c) == lo[c] && to_upper((char)c) == up[c] && swap_case((char)c) == sw[c]);  // 1 文字版
        STRB_CHECK(caesar_shift((char)c, 27) == caesar_shift(std::string(1, (char)c), 27)[0]);
    }
    for (int it = 0; it < 2000 * TEST_SCALE; ++it) {
        const int n = (int)(rng() % 20);
        std::string s = rnd_str(rng, n, 26);
        for (char& c : s)
            if (rng() % 5 == 0) c = (char)('A' + rng() % 26);
            else if (rng() % 9 == 0) c = ' ';
        const long long k = rng() % 2 ? (long long)(rng() % 60) - 30 : rnd_ll(rng);
        const std::string t = caesar_shift(s, k);
        for (int i = 0; i < n; ++i) {
            if ('a' <= s[i] && s[i] <= 'z') STRB_CHECK(t[i] == (char)('a' + (int)((s[i] - 'a' + ((k % 26) + 26) % 26) % 26)));
            else STRB_CHECK(t[i] == s[i]);
        }
        STRB_CHECK(caesar_shift(caesar_shift(s, k), -(k % 26)) == s);
        // CharMap：replace を何回か重ねる ⇔ 文字列に順に当てはめる
        CharMap cm;
        std::string bf = s;
        const int q = (int)(rng() % 8);
        for (int j = 0; j < q; ++j) {
            const char c = (char)('a' + rng() % 4), d = (char)('a' + rng() % 4);
            cm.replace(c, d);
            for (char& x : bf)
                if (x == c) x = d;
        }
        STRB_CHECK(cm.apply(s) == bf);
        CharMap g;
        g.set('a', 'z');
        g.replace('b', 'a');
        STRB_CHECK(cm.compose(g).apply(s) == g.apply(cm.apply(s)));
        ++cases;
    }
    // 整数の列の CharMap
    CharMap im(Alphabet::ints(-2, 3));
    im.replace(-2, 2), im.replace(2, 0);
    STRB_CHECK(im.apply(std::vector<int>{-2, -1, 0, 1, 2}) == std::vector<int>({0, -1, 0, 1, 0}) && im(7) == 7);
}

static void test_misc(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 5000 * TEST_SCALE; ++it) {
        const int n = (int)(rng() % 12), m = (int)(rng() % 12);
        const int sg = 1 + (int)(rng() % 3);
        const std::string a = rnd_str_shaped(rng, n, sg), b = rnd_str_shaped(rng, m, sg);
        std::string ra(a.rbegin(), a.rend());
        STRB_CHECK(reversed(a) == ra && min_with_reverse(a) == std::min(a, ra));
        std::string sa = a;
        std::sort(sa.begin(), sa.end());
        STRB_CHECK(sorted_chars(a) == sa);
        const long long k = rng() % 4;
        std::string rep;
        for (long long j = 0; j < k; ++j) rep += a;
        STRB_CHECK(repeat_string(a, k) == rep);
        int p = 0;
        while (p < n && p < m && a[p] == b[p]) ++p;
        STRB_CHECK(common_prefix_length(a, b) == p);
        int q = 0;
        while (q < n && q < m && a[n - 1 - q] == b[m - 1 - q]) ++q;
        STRB_CHECK(common_suffix_length(a, b) == q);
        if (n == m) {
            int h = 0;
            for (int i = 0; i < n; ++i) h += a[i] != b[i];
            STRB_CHECK(hamming_distance(a, b) == h);
        }
        // 部分列か：全部の部分集合で
        bool sub = false;
        if (n <= 10)
            for (int mask = 0; mask < (1 << m) && !sub; ++mask) {
                if (__builtin_popcount((unsigned)mask) != n) continue;
                std::string t;
                for (int i = 0; i < m; ++i)
                    if (mask >> i & 1) t += b[i];
                sub = t == a;
            }
        if (n <= 10) STRB_CHECK(is_subsequence(a, b) == sub);
        STRB_CHECK(is_palindrome(a) == (a == ra));
        // substr_of / StrSpan
        const int l = (int)(rng() % (n + 1)), r = l + (int)(rng() % (n - l + 1));
        STRB_CHECK(substr_of(a, StrSpan{l, r}) == a.substr(l, r - l) && StrSpan{l, r}.len() == r - l);
        // vector でも同じ
        const std::vector<int> va(a.begin(), a.end()), vb(b.begin(), b.end());
        STRB_CHECK(common_prefix_length(va, vb) == p && common_suffix_length(va, b) == q && is_subsequence(va, vb) == is_subsequence(a, b));
        STRB_CHECK(is_palindrome(va) == is_palindrome(a) && reversed(va) == std::vector<int>(ra.begin(), ra.end()));
        // split / join
        std::string sp = a;
        for (char& c : sp)
            if (c == 'a') c = ',';
        std::vector<std::string> parts, parts2;
        {
            std::string cur;
            for (char c : sp) {
                if (c == ',') parts.push_back(cur), cur.clear();
                else cur += c;
            }
            parts.push_back(cur);
            for (const auto& x : parts)
                if (!x.empty()) parts2.push_back(x);
        }
        STRB_CHECK(split(sp, ',') == parts && split(sp, ',', true) == parts2 && join(parts, ",") == sp);
        ++cases;
    }
    // 文字列リテラル・混在
    STRB_CHECK(is_subsequence("ace", std::string("abcde")) && !is_subsequence(std::string("aec"), "abcde"));
    STRB_CHECK(common_prefix_length("abc", "abd") == 2 && is_palindrome("racecar") && reversed("abc") == "cba");
    STRB_CHECK(count_chars("hello")[(int)('l' - 'a')] == 2 && run_length("aab").size() == 2);
    STRB_CHECK(join(std::vector<long long>{1, -2, 30}) == "1 -2 30" && join(std::vector<char>{'a', 'b'}, "") == "ab");
    STRB_CHECK(split("", ',') == std::vector<std::string>({""}) && split("", ',', true).empty());
    // 128 以上の char の比較（unsigned char として）
    const std::string hi = std::string(1, (char)200), lo2 = "a";
    STRB_CHECK(min_with_reverse(hi + lo2) == lo2 + hi && sorted_chars(hi + lo2) == lo2 + hi);
}

// 最終照合で足したもの
template <class V>
static bool bf_no_adjacent(const V& s) {
    for (std::size_t i = 1; i < s.size(); ++i)
        if (s[i] == s[i - 1]) return false;
    return true;
}
static void test_final_additions(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 3000 * TEST_SCALE; ++it) {
        const int n = (int)(rng() % 9), m = (int)(rng() % 9);
        const int sg = 1 + (int)(rng() % 3);
        std::string a = rnd_str_shaped(rng, n, sg), b = rnd_str_shaped(rng, m, sg);
        for (char& c : a)
            if (rng() % 5 == 0) c = '?';
        for (char& c : b)
            if (rng() % 5 == 0) c = '?';
        // ワイルドカードつきの共通接頭辞・接尾辞
        auto eq = [](char x, char y) { return x == y || x == '?' || y == '?'; };
        int p = 0, q = 0;
        while (p < n && p < m && eq(a[p], b[p])) ++p;
        while (q < n && q < m && eq(a[n - 1 - q], b[m - 1 - q])) ++q;
        STRB_CHECK(common_prefix_length(a, b, '?') == p && common_suffix_length(a, b, '?') == q);
        // 部分列の位置（貪欲）
        const auto pos = subsequence_positions(a, b);
        STRB_CHECK(pos.empty() ? !is_subsequence(a, b) || a.empty() : (pos.size() == a.size() && is_subsequence(a, b)));
        for (std::size_t k = 0; k < pos.size(); ++k) STRB_CHECK(b[pos[k]] == a[k] && (k == 0 || pos[k - 1] < pos[k]));
        if (!pos.empty()) {  // 貪欲：それぞれ前の位置の次にある最初の一致
            int from = 0;
            for (std::size_t k = 0; k < pos.size(); ++k) {
                int j = from;
                while (b[j] != a[k]) ++j;
                STRB_CHECK(pos[k] == j);
                from = j + 1;
            }
        }
        // 隣が同じにならない並べ替え
        const std::string base = rnd_str_shaped(rng, n, sg);
        std::string mn;
        std::string perm = base;
        std::sort(perm.begin(), perm.end());
        do {
            if (bf_no_adjacent(perm) && (mn.empty() || perm < mn)) mn = perm;
        } while (std::next_permutation(perm.begin(), perm.end()));
        const std::string any = rearrange_no_adjacent_equal(base), lex = lex_smallest_no_adjacent_equal(base);
        STRB_CHECK(lex == mn && (mn.empty() ? any.empty() : (bf_no_adjacent(any) && sorted_chars(any) == sorted_chars(base))));
        // 隣どうしの条件で区切る
        const auto runs = split_runs_by(base, [](char x, char y) { return x != y; });
        int cur = 0;
        for (const StrSpan& sp : runs) {
            STRB_CHECK(sp.l == cur && sp.len() > 0 && bf_no_adjacent(substr_of(base, sp)));
            if (sp.r < n) STRB_CHECK(base[sp.r] == base[sp.r - 1]);
            cur = sp.r;
        }
        STRB_CHECK(cur == n);
        // 区切り（文字列）・前後の文字を除く
        std::string t = a;
        for (char& c : t)
            if (c == '?') c = ' ';
        std::string tt;
        for (char c : t) tt += c == ' ' ? std::string(", ") : std::string(1, c);
        STRB_CHECK(split(tt, ", ") == split(t, ' ') && split(tt, std::string(", "), true) == split(t, ' ', true));
        std::size_t L = 0, Rr = t.size();
        while (L < Rr && t[L] == ' ') ++L;
        while (Rr > L && t[Rr - 1] == ' ') --Rr;
        std::size_t R2 = t.size();
        while (R2 > 0 && t[R2 - 1] == ' ') --R2;
        STRB_CHECK(strip(t) == t.substr(L, Rr - L) && strip_left(t, ' ') == t.substr(L) && strip_right(t) == t.substr(0, R2));
        ++cases;
    }
    // 数の列でも
    for (int it = 0; it < 500 * TEST_SCALE; ++it) {
        std::vector<long long> v(rng() % 8);
        for (auto& x : v) x = rng() % 3 == 0 ? rnd_ll(rng) % 2 : (long long)(rng() % 3) - 1;
        std::vector<long long> mn, perm = v;
        std::sort(perm.begin(), perm.end());
        bool found = false;
        do {
            if (bf_no_adjacent(perm) && (!found || perm < mn)) mn = perm, found = true;
        } while (std::next_permutation(perm.begin(), perm.end()));
        const auto lex = lex_smallest_no_adjacent_equal(v), any = rearrange_no_adjacent_equal(v);
        STRB_CHECK(found ? (lex == mn && bf_no_adjacent(any) && sorted_chars(any) == sorted_chars(v)) : (v.empty() || (lex.empty() && any.empty())));
        ++cases;
    }
    // 大きな入力：全部同じ数 + 1
    {
        std::string big(100000, 'a');
        big += std::string(99999, 'b');
        const std::string r = lex_smallest_no_adjacent_equal(big);
        STRB_CHECK(r.size() == big.size() && r[0] == 'a' && bf_no_adjacent(r) && rearrange_no_adjacent_equal(big + "a").empty());
    }
    // 型・多重定義
    STRB_CHECK(from_ids(std::vector<long long>{7, 4, 11, 11, 14}) == "hello" && from_ids(std::vector<int>{0, 1}, Alphabet::upper()) == "AB");
    STRB_CHECK(join(std::set<int>{3, 1, 2}) == "1 2 3" && join(std::vector<std::string>{"a", "b"}, '\n') == "a\nb");
    STRB_CHECK(Alphabet::lower().with("@#").id('#') == 27 && count_chars("a@b", Alphabet::lower().with("@"))[26] == 1);
    STRB_CHECK(is_subsequence('b', "abc") && common_prefix_length('a', "ab") == 1);  // 1 文字も列として
    STRB_CHECK(split("a, b, , c", ", ") == std::vector<std::string>({"a", "b", "", "c"}) && split("", ", ").size() == 1);
    STRB_CHECK(strip("**a*b**", '*') == "a*b" && strip_right("ab  ") == "ab" && strip("   ").empty() && strip_right("ab\r\n") == "ab");
    STRB_CHECK(join(std::vector<double>{1234567.5, 0.1}) == "1234567.500000000000000 0.100000000000000" && join(std::vector<const char*>{"x", "y"}) == "x y");
    const auto cc = count_chars(std::string(100000, 'a'));
    STRB_CHECK(cc[0] * (cc[0] - 1) / 2 == 4999950000LL);  // long long なので積が溢れない
}

static void random_test() {
    std::mt19937_64 rng(20261012);
    int cases = 0;
    test_alphabet(rng, cases);
    test_compress(rng, cases);
    test_char_queries(rng, cases);
    test_run_length(rng, cases);
    test_transform(rng, cases);
    test_misc(rng, cases);
    test_final_additions(rng, cases);
    std::cout << "[str_base] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[str_base] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    auto ms = [](std::chrono::steady_clock::time_point a) {
        return (long long)std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count();
    };
    const int n = 1000000, Q = 10000000;
    const std::string s = rnd_str(rng, n, 26);
    std::vector<int> ls(Q), rs(Q);
    for (int i = 0; i < Q; ++i) {
        ls[i] = (int)(rng() % (n + 1)), rs[i] = (int)(rng() % (n + 1));
        if (ls[i] > rs[i]) std::swap(ls[i], rs[i]);
    }
    auto t0 = std::chrono::steady_clock::now();
    const CharPrefix cp(s);
    const long long t_cp_build = ms(t0);
    t0 = std::chrono::steady_clock::now();
    long long s1 = 0;
    for (int i = 0; i < Q; ++i) s1 += cp.count(ls[i], rs[i], (char)('a' + (i % 26)));
    const long long t_cp_q = ms(t0);
    t0 = std::chrono::steady_clock::now();
    const CharPositions cps(s);
    const long long t_ps_build = ms(t0);
    t0 = std::chrono::steady_clock::now();
    long long s2 = 0;
    for (int i = 0; i < Q / 10; ++i) s2 += cps.next(ls[i], (char)('a' + (i % 26))) + cps.prev(rs[i], (char)('a' + (i % 26)));
    const long long t_ps_q = ms(t0);
    const int n2 = 200000;
    const CharPositions cps2(s.substr(0, n2));
    t0 = std::chrono::steady_clock::now();
    for (int i = 0; i < Q / 10; ++i) s2 += cps2.next(ls[i] % (n2 + 1), (char)('a' + (i % 26))) + cps2.prev(rs[i] % (n2 + 1), (char)('a' + (i % 26)));
    const long long t_ps_q2 = ms(t0);
    t0 = std::chrono::steady_clock::now();
    for (int i = 0; i < Q / 10; ++i) s2 += cp.next(ls[i], (char)('a' + (i % 26))) + cp.prev(rs[i], (char)('a' + (i % 26)));
    const long long t_cp_np = ms(t0);
    const std::string big = rnd_str_shaped(rng, 10000000, 3);
    t0 = std::chrono::steady_clock::now();
    const auto runs = run_length(big);
    const long long t_rle = ms(t0);
    t0 = std::chrono::steady_clock::now();
    const auto ids = to_ids(big);
    const long long t_ids = ms(t0);
    CharMap cm;
    cm.replace('a', 'b'), cm.replace('c', 'a');
    t0 = std::chrono::steady_clock::now();
    const std::string mapped = cm.apply(big);
    const long long t_map = ms(t0);
    std::vector<long long> a(1000000);
    for (auto& x : a) x = rnd_ll(rng);
    t0 = std::chrono::steady_clock::now();
    const auto comp = compress_values(a);
    const long long t_comp = ms(t0);
    std::cout << "[str_base] benchmark: CharPrefix（N = 10^6, σ = 26）構築 " << t_cp_build << " ms・count × 10^7 " << t_cp_q
              << " ms・next + prev × 10^6 " << t_cp_np << " ms\n"
              << "[str_base] benchmark: CharPositions（N = 10^6）構築 " << t_ps_build << " ms・next + prev × 10^6 " << t_ps_q << " ms（N = 2×10^5 なら " << t_ps_q2 << " ms）\n"
              << "[str_base] benchmark: run_length（N = 10^7）" << t_rle << " ms, to_ids（N = 10^7）" << t_ids << " ms, CharMap::apply（N = 10^7）"
              << t_map << " ms, compress_values（N = 10^6）" << t_comp << " ms  (check "
              << (s1 + s2 + (long long)runs.size() + ids[5] + mapped[7] + comp.size()) % 7 << ")\n";
}

}  // namespace strb_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace strb_demo {

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
#define SHOW(expr, comment) strb_demo::show(#expr, (expr), comment)
#define RUN(...)                   \
    strb_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static std::string runs_to_s(const std::vector<RunBlock<char>>& rs) {
    std::string r = "[";
    for (std::size_t i = 0; i < rs.size(); ++i)
        r += (i ? ", " : "") + std::string("{'") + rs[i].value + "', " + std::to_string(rs[i].len) + ", " + std::to_string(rs[i].start) + "}";
    return r + "]";
}

static void demo() {
    banner("文字列基盤（str_base） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  Alphabet::lower / upper / digits / letters / alnum / ascii / bytes / of / ints / with   文字の種類（番号の順 = 辞書順）\n"
              << "  to_ids / from_ids / compress_values             文字列 ↔ 番号の列 / 数の列を 0〜K-1 に\n"
              << "  count_chars / CharPrefix / CharPositions        文字の数 / 区間の文字の数 / 次・前の出現・k 番目（CharPrefix は番号の版 count_id なども）\n"
              << "  run_length / run_length_decode                  ランレングス圧縮 / 元に戻す\n"
              << "  to_lower / to_upper / swap_case / caesar_shift / CharMap   変換（1 文字の版も）/ 文字の置き換えを重ねる\n"
              << "  reversed / min_with_reverse / sorted_chars / repeat_string  反転 / 向きを区別しない代表 / 並べ替え / 繰り返し\n"
              << "  common_prefix_length / common_suffix_length / hamming_distance / is_subsequence / subsequence_positions / is_palindrome  比較\n"
              << "  split / join / strip / strip_left / strip_right / StrSpan / substr_of   区切る / つなげる（答えの出力にも）/ 両端を除く / 部分文字列の位置\n"
              << "  split_runs_by / rearrange_no_adjacent_equal / lex_smallest_no_adjacent_equal   隣どうしの条件で区切る / 隣が同じにならない並べ替え\n"
              << "  共通：範囲は [l, r)（L 文字目〜R 文字目は (L - 1, R)）。位置がないときは -1。string も vector<int> も渡せる\n";

    api("【1】 Alphabet（文字の種類）", "構造体",
        "文字・値に 0 からの番号を付ける。ほかの文字列テンプレート（Trie・NextTable など）にも渡す。範囲外の文字はメッセージを出して止まる", "O(1)");
    RUN(const Alphabet L = Alphabet::lower());
    SHOW(L.size(), "a〜z の 26 種類（既定）");
    SHOW(L.id('c'), "'c' の番号");
    SHOW(L.ch(25), "番号 25 の文字");
    SHOW(L.contains('A'), "大文字は入っていない（id('A') は止まる）");
    SHOW(Alphabet::letters().id('a'), "A〜Z が 0〜25、a〜z が 26〜51（ASCII 順）");
    RUN(const Alphabet B = Alphabet::of("()[]"));
    SHOW(B.id(']'), "並べた順に番号（並べた順が辞書順になる）");
    SHOW(to_ids("ba", Alphabet::of("bac")) < to_ids("ab", Alphabet::of("bac")), "独自の順 b < a < c で \"ba\" < \"ab\"（ABC219 C）");
    RUN(const Alphabet K = Alphabet::ints(-5, 5));
    SHOW(K.id(-5), "整数 [-5, 5) の番号は 0〜9");
    SHOW(K.value(9), "番号 → 値");
    SHOW(to_ids("cab"), "文字列 → 番号の列");
    SHOW(from_ids({7, 8}), "番号の列 → 文字列（vector<long long> でもよい）");
    SHOW(Alphabet::lower().with("@").id('@'), "英小文字のあとに '@' を足した種類（27 種類）");
    note("char は 0〜255 として扱う（日本語などの 128 以上の文字は Alphabet::bytes() で）。");
    note("既定は英小文字なので、大文字・数字・記号が混ざると「文字の種類に入っていない」で止まる → upper() / digits() / bytes() / with(\"…\") を渡す。");

    api("【2】 compress_values(a)", "CompressedValues{ids, values}",
        "数の列を大小の順を保って 0〜K-1 にする。.alphabet() が Alphabet::ints(K) を返すので、値が大きい列も Trie などに渡せる", "O(N log N)");
    RUN(const auto c = compress_values(std::vector<long long>{1000000000000LL, -7, 1000000000000LL, 3}));
    SHOW(c.ids, "各要素の番号");
    SHOW(c.values, "番号 → 元の値（昇順）");
    SHOW(c.alphabet().size(), "K = 3");

    api("【3】 count_chars(s, A) / CharPrefix cp(s, A)", "vector<long long> / 構造体（数は long long）",
        "文字の個数。CharPrefix は区間 [l, r) の文字の数と、i 以上で最初の c・i より前で最後の c を O(1) で（前処理 O(Nσ)、メモリ (N+1)σ 個の int）", "O(N) / O(1)");
    RUN(const std::string S = "abracadabra");
    SHOW(count_chars(S), "a〜z の個数");
    RUN(const CharPrefix cp(S));
    SHOW(cp.count(0, 4, 'a'), "\"abra\" の 'a' の数");
    SHOW(cp.has(5, 7, 'r'), "\"da\" に 'r' はない");
    SHOW(cp.distinct(0, 11), "文字の種類の数（a, b, c, d, r）");
    SHOW(cp.next(4, 'b'), "4 以上で最初の 'b'");
    SHOW(cp.prev(4, 'b'), "4 より前で最後の 'b'");
    note("問題の「L 文字目から R 文字目」は cp.count(L - 1, R, c)。");
    SHOW(cp.count_id(0, 11, 0), "番号で（rep(k, 26) で回すときは count_id / total_id / next_id / prev_id / kth_id。cp.total(k) は値 k の文字になる）");
    note("記号などを全部数えたいときは count_chars(S, Alphabet::bytes())[c]（番号 = 文字コード）。");
    note("N σ が大きすぎる（メモリが足りない）ときは、次の CharPositions（O(log N)・メモリ O(N + σ)）。");

    api("【4】 CharPositions cps(s, A)", "構造体",
        "文字ごとの出現位置（CharPrefix のメモリを抑えた版）。i 以上で最初の c・i より前で最後の c・c の k 番目。メモリ O(N + σ) なので σ が大きくてもよい", "O(log N) / O(1)");
    RUN(const CharPositions cps(S));
    SHOW(cps.next(4, 'a'), "4 以上で最初の 'a'");
    SHOW(cps.prev(4, 'a'), "4 より前で最後の 'a'");
    SHOW(cps.next(11, 'a'), "ない → -1");
    SHOW(cps.kth('a', 2), "'a' の 2 番目（0 から）の位置");
    SHOW(cps.positions('r'), "'r' の位置の一覧");
    SHOW(cps.count(0, 8, 'a'), "[0, 8) の 'a' の数");
    SHOW(cps.occurrence_index(7), "S[7] = 'a' は 'a' の何番目の出現か");

    api("【5】 run_length(s) / run_length_decode(runs)", "vector<RunBlock<T>{value, len, start}> / 元の列",
        "連続する同じ値をまとめる。len は long long。vector<long long> でもよい", "O(N)");
    std::cout << "  " << pad("run_length(\"aaabccdddd\")", 34) << " -> " << runs_to_s(run_length("aaabccdddd")) << "\n";
    SHOW(run_length_decode(run_length("aaabccdddd")), "元に戻す");
    {
        // ABC259 C: S の「同じ文字が 2 つ以上続くところ」にその文字を足していき、T にできるか
        const std::string S2 = "abbaac", T2 = "abbbbaaac";
        const auto a = run_length(S2), b = run_length(T2);
        bool ok = a.size() == b.size();
        for (std::size_t i = 0; ok && i < a.size(); ++i)
            ok = a[i].value == b[i].value && (a[i].len == b[i].len || (a[i].len >= 2 && a[i].len <= b[i].len));
        note("ABC259 C（XX to XXX）: 両方を run_length して、文字が同じ・長さが同じか「S が 2 以上で T 以下」か");
        std::cout << "    S = abbaac, T = abbbbaaac → " << (ok ? "Yes" : "No") << "\n";
    }
    {
        // 典型90 084: 'o' と 'x' の両方を含む部分文字列の数 = 全部 − 1 種類だけのもの（run ごとに len(len+1)/2）
        const std::string S3 = "ooxo";
        long long ans = (long long)S3.size() * ((long long)S3.size() + 1) / 2;
        for (const auto& r : run_length(S3)) ans -= r.len * (r.len + 1) / 2;
        note("典型90 084（o と x の両方を含む区間の数）: 全区間 − Σ len(len+1)/2");
        std::cout << "    S = ooxo → " << ans << "\n";
    }
    {
        // ABC329 C: 1 種類の文字だけからなる部分文字列（文字列として）の数 = 文字ごとの最長の run の和
        std::vector<long long> best(26, 0);
        for (const auto& r : run_length("aaabaa")) best[r.value - 'a'] = std::max(best[r.value - 'a'], r.len);
        note("ABC329 C（Count xxx）: 文字ごとに最長の run の長さを足す");
        std::cout << "    S = aaabaa → " << std::accumulate(best.begin(), best.end(), 0LL) << "\n";
    }

    api("【6】 to_lower / to_upper / swap_case / caesar_shift(s, k, A)", "string",
        "英字の大文字・小文字の変換。caesar_shift は A の文字を番号で k ずらす（負も。A にない文字はそのまま）", "O(N)");
    SHOW(to_upper("AtCoder"), "");
    SHOW(swap_case("AtCoder"), "");
    SHOW(caesar_shift("xyz abc", 3), "a〜z を 3 ずらす（空白はそのまま）");
    SHOW(caesar_shift("HAL", 1, Alphabet::upper()), "大文字を 1 ずらす");
    SHOW(to_upper('q'), "1 文字の版（to_lower・swap_case・caesar_shift(c, k, A) も）");
    {
        // ABC232 B: S をシーザー暗号でずらして T にできるか
        const std::string S4 = "abc", T4 = "ijk";
        bool ok = false;
        for (int k = 0; k < 26; ++k) ok |= caesar_shift(S4, k) == T4;
        note("ABC232 B（Caesar Cipher）: k = 0〜25 で caesar_shift(S, k) == T を確かめる");
        std::cout << "    S = abc, T = ijk → " << (ok ? "Yes" : "No") << "\n";
    }

    api("【7】 CharMap cm(A)：replace(c, d) / set(c, d) / cm(c) / apply(s) / compose(g)", "構造体",
        "「今 c になっている文字を全部 d に」を何度も重ねる（毎回文字列を書き換えると O(NQ)。表だけ持てば O(σQ + N)）", "replace O(σ)・apply O(N)");
    {
        // ABC342 C: 「c を全部 d に」を Q 回
        RUN(CharMap cm);
        RUN(cm.replace('r', 'a'));  // atcoder → atcodea
        RUN(cm.replace('t', 'e'));  // → aecodea
        RUN(cm.replace('d', 'v'));  // → aecovea
        RUN(cm.replace('a', 'r'));  // 元の 'a' も元の 'r' も 'r' に → recover
        SHOW(cm.apply("atcoder"), "ABC342 C の入力例 1");
        SHOW(cm('r'), "元の 'r' の行き先（a を経て r）");
    }

    api("【8】 reversed / min_with_reverse / sorted_chars / repeat_string", "string（vector なら vector）",
        "反転 / s と反転の小さい方（向きを区別しないときの代表）/ 並べ替え（並べ替えて同じかの比較用）/ k 回つなげる", "O(N)");
    SHOW(reversed("abc"), "");
    SHOW(min_with_reverse("cba"), "\"abc\" と同じものとみなす");
    {
        // ABC310 C: 反転したものを同じとみなして、種類を数える
        const std::vector<std::string> sticks = {"a", "abc", "de", "cba", "de", "abc"};
        std::set<std::string> st;
        for (const auto& x : sticks) st.insert(min_with_reverse(x));
        note("ABC310 C（Reversible）: set に min_with_reverse(S_i) を入れて数える");
        std::cout << "    {a, abc, de, cba, de, abc} → " << st.size() << "\n";
    }
    SHOW(sorted_chars("listen") == sorted_chars("silent"), "並べ替えて同じ（アナグラム）");
    SHOW(repeat_string("ab", 3), "");

    api("【9】 common_prefix_length / common_suffix_length / hamming_distance / is_subsequence(t, s) / is_palindrome(s)", "int / bool",
        "共通の接頭辞・接尾辞の長さ、同じ位置で違う文字の数、t が s の部分列（飛び飛び）か、全体が回文か", "O(N)");
    SHOW(common_prefix_length("atcoder", "atcorder"), "");
    SHOW(common_suffix_length("atcoder", "beginner"), "\"er\"");
    SHOW(hamming_distance("karolin", "kathrin"), "");
    SHOW(is_subsequence("acd", "abcde"), "a・c・d を順に拾える");
    SHOW(is_subsequence("dc", "abcde"), "順番が逆");
    SHOW(is_palindrome("level"), "");
    SHOW(is_palindrome(std::vector<int>{1, 2, 1}), "数の列でも");
    SHOW(common_prefix_length("ab?de", "a?cdx", '?'), "'?' はどちらにあっても一致（ABC287 D）");
    SHOW(subsequence_positions("ace", "abcde"), "t を s に左から貪欲に当てはめた位置（部分列でなければ空。ABC352 B）");
    note("is_subsequence(t, s) は「t が s の部分列か」。短い方（拾う方）を先に書く。");
    note("同じ文字列に何度も部分列か聞くなら 06_subsequence の NextTable、区間の回文判定は 04_palindrome。");

    api("【10】 split(s, sep, skip_empty) / join(vs, sep)", "vector<string> / string",
        "sep で区切る（空の部分も残す。skip_empty = true で除く）/ 区切りを挟んでつなげる（数の列も。答えの出力に）", "O(N)");
    SHOW(split("2026/10/12", '/'), "");
    SHOW(split("a,,b", ','), "空の部分も残る");
    SHOW(split("a,,b", ',', true), "空の部分を除く");
    SHOW(join(std::vector<std::string>{"x", "y", "z"}, "-"), "");
    SHOW(join(std::vector<int>{3, 1, 4}), "数の列を空白区切りで（cout << join(ans) << '\\n';）");
    SHOW(split("a, b, c", ", "), "文字列の区切り");
    SHOW(join(std::set<int>{5, 2}, '\n') == "2\n5", "set なども・区切りを 1 文字で");
    SHOW(strip_right("ab***", '*'), "末尾の '*' を除く（strip_left・strip も。既定は空白）");

    api("【11】 StrSpan{l, r} / substr_of(s, sp)", "構造体 / string",
        "部分文字列の位置 [l, r)（len() で長さ）。最長の回文・k 番目の部分文字列などはこの形で返す", "O(1) / O(長さ)");
    RUN(const StrSpan sp{2, 5});
    SHOW(sp.len(), "");
    SHOW(substr_of("abcdefg", sp), "[2, 5) の部分");

    api("【12】 split_runs_by(s, same) / rearrange_no_adjacent_equal(s) / lex_smallest_no_adjacent_equal(s)", "vector<StrSpan> / string（vector なら vector）",
        "same(s[i-1], s[i]) が成り立つ間をまとめた区間 / 隣どうしが違う並べ替え（できなければ空）/ そのうち辞書順最小", "O(N) / O(N log σ)");
    SHOW(split_runs_by("abbcab", [](char x, char y) { return x != y; }), "隣が違う区間に分ける（ABC456 C）");
    SHOW(rearrange_no_adjacent_equal("aaabb"), "ABC459 D");
    SHOW(lex_smallest_no_adjacent_equal("aabbc"), "");
    SHOW(rearrange_no_adjacent_equal("aaab").empty(), "一番多い文字が (N+1)/2 個を超えるとできない");

    head("【型】 巨大な文字列の k 文字目（ABC380 D Strange Mirroring の形）");
    note("「S の大文字と小文字を入れ替えた T を後ろにつなげる」を 10^100 回 → K 文字目（1 始まり）は？");
    note("k = K - 1 を |S| で割った商 q の 1 のビットの数が奇数なら、S[k mod |S|] の大文字・小文字が入れ替わっている。");
    {
        const std::string S5 = "aB";
        std::string out;
        for (long long K1 = 1; K1 <= 16; ++K1) {
            const long long k = K1 - 1, q = k / (long long)S5.size();
            const char ch = S5[k % (long long)S5.size()];
            out += (__builtin_popcountll((unsigned long long)q) & 1) ? swap_case(std::string(1, ch))[0] : ch;
        }
        std::cout << "    S = aB の 1〜16 文字目 → " << out << "（実際に倍々にした " << [&] {
            std::string t = S5;
            while (t.size() < 16) t += swap_case(t);
            return t.substr(0, 16);
        }() << " と同じ）\n";
    }
    head("【型】 1 文字が 2 文字に増える変換を t 回（ABC242 D ABC Transform の形）");
    note("A → BC、B → CA、C → AB を t 回。k 文字目（0 から）の親は k / 2 文字目で、子は親を (k % 2 + 1) ずらした文字。");
    note("k が 0 になるか t が 0 になるまで親をたどり、ずれの合計で caesar_shift(c, ずれ, Alphabet::of(\"ABC\"))。");
    {
        auto at = [](const std::string& S0, long long t, long long k) {
            long long sh = 0;
            while (t > 0 && k > 0) sh += k % 2 + 1, k /= 2, --t;
            return caesar_shift(S0[(std::size_t)k], sh + t, Alphabet::of("ABC"));
        };
        std::cout << "    S = ABC：(t, k) = (0, 0), (1, 0), (1, 2), (1, 5) → " << at("ABC", 0, 0) << at("ABC", 1, 0) << at("ABC", 1, 2) << at("ABC", 1, 5) << "（入力例 1：A B C B）\n";
    }
    head("【型】 S_n = S_{n-1} + S_{n-2} の区間の文字の数（ABC450 E Fibonacci String の形）");
    note("長さが 10^18 を超えるまで「長さ」と「文字ごとの数」だけを作り、先頭 r 文字の数を n の大きい方から分けて求める。");
    note("土台の X・Y の接頭辞は CharPrefix の count。区間 [L, R] は pre(R) − pre(L − 1)。");
}

}  // namespace strb_demo

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
    strb_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    strb_test::benchmark();    // 速度計測
    strb_demo::demo();         // 使い方デモ
    return 0;
}
