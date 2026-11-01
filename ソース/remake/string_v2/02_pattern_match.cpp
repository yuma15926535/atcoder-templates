// =============================================================================
//  パターン検索（Z・KMP・出現の列挙と数・置き換え・周期・境界・回転・重なり・KMP オートマトン）   --- 文字列第二世代
// -----------------------------------------------------------------------------
//  分類       : 02 一致・出現
//  依存       : 00_str_base（同梱済み）
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）
//  文字コード : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  text・pat は string でも vector<int>・vector<long long> でもよい（混ぜてもよい。文字列リテラル "abc" も）
//  位置は 0 始まり。「ない」ときは -1。ACL の z_algorithm（z[0] = |s|）と同じ約束で、名前は z_array
//
//  機能一覧                                                                       計算量
//    z_array(s)                        z[i] = s と s[i..) の最長共通接頭辞の長さ（z[0] = |s|）       O(N)
//    prefix_function(s)                π[i] = s[0, i] の最長の真の境界（接頭辞 = 接尾辞）の長さ       O(N)
//    match_lengths(text, pat)          各 i で text[i..) と pat の最長共通接頭辞の長さ                O(N + M)
//    find_all(text, pat)               pat が現れる開始位置（重なりも。昇順）                         O(N + M)
//    count_occurrences(text, pat)      現れる回数（重なりも）                                         O(N + M)
//    find_first(text, pat, from) / find_last(text, pat)   from 以上で最初 / 最後の位置（なければ -1）   O(N + M)
//    count_non_overlapping(text, pat)  左から重ならないように取った数                                 O(N + M)
//    replace_all(text, from, to)       from を左から重ならないように全部 to に                         O(N + M)
//    all_borders(s)                    0 < k < |s| で接頭辞 = 接尾辞となる長さ k（大きい順）            O(N)
//    min_period(s) / all_periods(s)    最小の周期 / すべての周期（|s| も含む。小さい順）                O(N)
//    repetition_root(s)                s = t^k となる最短の t の長さ                                  O(N)
//    prefix_min_periods(s)             各 i で s[0, i) の最小の周期（長さ |s| + 1）                     O(N)
//    rotation_offsets(a, b) / is_rotation(a, b)   a を左に k 回まわすと b になる k の一覧 / なるか     O(N)
//    find_all_cyclic(text, pat)        text を無限に繰り返した中で pat が始まる位置（0 ≤ i < |text|）   O(N + M)
//    max_overlap(a, b)                 a の接尾辞 = b の接頭辞となる最大の長さ                        O(N + M)
//    KmpAutomaton KA(pat, A)           状態 0〜|pat|（先頭何文字まで一致しているか）の遷移表          構築 O(Mσ)
//      KA.next(state, c) / accepting(state) / size()   遷移 / pat が現れたところか / 状態の数       O(1)
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    pm_test::random_test()  … 愚直な比較（全部の位置で比べる）と全機能を突き合わせるランダムテスト
//    pm_test::benchmark()    … 速度計測
//    pm_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
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
#ifndef GG_PATTERN_MATCH
#define GG_PATTERN_MATCH
// ---- パターン検索（文字列第二世代）  Z・KMP・出現・周期・境界・回転・重なり・置き換えの繰り返し・区間の出現・KMP オートマトン ----
// text・pat は string でも vector<整数> でもよい（混ぜてもよい。値で比べる）。文字列リテラルも渡せる
namespace str_internal {
// 列（string・vector<整数>）をそのまま読む（コピーしない）。要素は code_of で比べる
template <class S>
std::vector<int> z_codes(const S& s) {
    const int n = (int)s.size();
    std::vector<int> z(n);
    if (n == 0) return z;
    z[0] = n;
    for (int i = 1, l = 0, r = 0; i < n; ++i) {
        int k = i < r ? std::min(r - i, z[i - l]) : 0;
        while (i + k < n && s[k] == s[i + k]) ++k;
        z[i] = k;
        if (i + k > r) l = i, r = i + k;
    }
    return z;
}
template <class S>
std::vector<int> pi_codes(const S& s) {
    const int n = (int)s.size();
    std::vector<int> pi(n);
    for (int i = 1; i < n; ++i) {
        int k = pi[i - 1];
        while (k > 0 && !(s[i] == s[k])) k = pi[k - 1];
        pi[i] = k + (s[i] == s[k]);
    }
    return pi;
}
// text[from..) の中の pat の出現開始位置（text での位置。重なりも）。pat が空なら from〜|text|
template <class T, class P>
std::vector<int> kmp_find(const T& t, const P& p, int from = 0) {
    const int n = (int)t.size(), m = (int)p.size();
    std::vector<int> r;
    if (m == 0) {
        for (int i = from; i <= n; ++i) r.push_back(i);
        return r;
    }
    const std::vector<int> pi = pi_codes(p);
    for (int i = from, k = 0; i < n; ++i) {
        const long long c = code_of(t[i]);
        while (k > 0 && (k == m || c != code_of(p[k]))) k = pi[k - 1];
        if (c == code_of(p[k])) ++k;
        if (k == m) r.push_back(i - m + 1);
    }
    return r;
}
}  // namespace str_internal

template <class Seq>
std::vector<int> z_array(const Seq& s) { return str_internal::z_codes(str_internal::seq(s)); }  // z[i] = s と s[i..) の最長共通接頭辞（z[0] = |s|）  O(N)
template <class Seq>
std::vector<int> prefix_function(const Seq& s) { return str_internal::pi_codes(str_internal::seq(s)); }  // π[i] = s[0, i] の最長の真の境界の長さ  O(N)
template <class T, class P>
std::vector<int> match_lengths(const T& text, const P& pat) {  // 各 i で text[i..) と pat の最長共通接頭辞の長さ（長さ |text|）  O(N + M)
    const auto& t = str_internal::seq(text);
    const auto& p = str_internal::seq(pat);
    const int n = (int)t.size(), m = (int)p.size();
    const std::vector<int> zp = str_internal::z_codes(p);
    std::vector<int> res(n);
    for (int i = 0, l = 0, r = 0; i < n; ++i) {  // [l, r) は text[l, r) == pat[0, r - l) となる右端が最大の区間
        int k = 0;
        if (i < r) k = std::min(r - i, zp[i - l]);
        if (i >= r || zp[i - l] >= r - i) {
            while (k < m && i + k < n && str_internal::code_of(t[i + k]) == str_internal::code_of(p[k])) ++k;
        }
        res[i] = k;
        if (i + k > r) l = i, r = i + k;
    }
    return res;
}
template <class T, class P>
std::vector<int> find_all(const T& text, const P& pat) {  // pat が現れる開始位置（重なりも。昇順。pat が空なら 0〜|text|）  O(N + M)
    return str_internal::kmp_find(str_internal::seq(text), str_internal::seq(pat));
}
template <class T, class P>
long long count_occurrences(const T& text, const P& pat) { return (long long)find_all(text, pat).size(); }  // 重なりも数える
template <class T, class P>
int find_first(const T& text, const P& pat, int from = 0) {  // from 以上で最初に現れる位置（なければ -1）  O(N + M)
    const auto& t = str_internal::seq(text);
    if (from < 0 || from > (int)t.size()) str_internal::fail("find_first", "from は 0〜|text| で指定する");
    const auto r = str_internal::kmp_find(t, str_internal::seq(pat), from);
    return r.empty() ? -1 : r[0];
}
template <class T, class P>
int find_last(const T& text, const P& pat) {  // 最後に現れる位置（なければ -1）
    const auto r = find_all(text, pat);
    return r.empty() ? -1 : r.back();
}
template <class T, class P>
long long count_non_overlapping(const T& text, const P& pat) {  // 左から重ならないように取った数（pat は空でない）
    const auto& p = str_internal::seq(pat);
    if (p.empty()) str_internal::fail("count_non_overlapping", "pat が空");
    long long c = 0, last = -1;
    for (int i : str_internal::kmp_find(str_internal::seq(text), p))
        if (i >= last) ++c, last = i + (long long)p.size();
    return c;
}
template <class T, class F, class U>
str_internal::seq_t<T> replace_all(const T& text0, const F& from, const U& to0) {  // from を左から重ならないように全部 to に（from は空でない）  O(N + M)
    const auto& text = str_internal::seq(text0);
    const auto& to = str_internal::seq(to0);
    const auto& f = str_internal::seq(from);
    if (f.empty()) str_internal::fail("replace_all", "from が空");
    str_internal::seq_t<T> r;
    std::size_t pos = 0;
    for (int i : str_internal::kmp_find(text, f)) {
        if ((std::size_t)i < pos) continue;
        r.insert(r.end(), text.begin() + pos, text.begin() + i);
        for (std::size_t j = 0; j < to.size(); ++j) r.push_back((typename str_internal::seq_t<T>::value_type)str_internal::code_of(to[j]));
        pos = (std::size_t)i + f.size();
    }
    r.insert(r.end(), text.begin() + pos, text.end());
    return r;
}

// ---- 周期・境界 ----
template <class Seq>
std::vector<int> all_borders(const Seq& s) {  // 0 < k < |s| で s[0, k) == s[|s|-k, |s|) となる k（大きい順）  O(N)
    const auto pi = prefix_function(s);
    std::vector<int> r;
    if (pi.empty()) return r;
    for (int k = pi.back(); k > 0; k = pi[k - 1]) r.push_back(k);
    return r;
}
template <class Seq>
int min_period(const Seq& s) {  // 最小の周期 p（s[i] == s[i + p]。割り切れなくてよい。空なら 0）  O(N)
    const auto pi = prefix_function(s);
    return pi.empty() ? 0 : (int)pi.size() - pi.back();
}
template <class Seq>
std::vector<int> all_periods(const Seq& s) {  // 1〜|s| のすべての周期（小さい順。|s| 自身も含む）  O(N)
    const int n = (int)str_internal::seq(s).size();
    std::vector<int> r;
    for (int k : all_borders(s)) r.push_back(n - k);
    if (n > 0) r.push_back(n);
    return r;
}
template <class Seq>
int repetition_root(const Seq& s) {  // s = t^k（t を k 回）となる最短の t の長さ（空なら 0）  O(N)
    const int n = (int)str_internal::seq(s).size(), p = min_period(s);
    return (n > 0 && n % p == 0) ? p : n;
}
template <class Seq>
std::vector<int> prefix_min_periods(const Seq& s) {  // r[i] = s[0, i) の最小の周期（長さ |s| + 1。r[0] = 0）  O(N)
    const auto pi = prefix_function(s);
    std::vector<int> r(pi.size() + 1, 0);
    for (std::size_t i = 0; i < pi.size(); ++i) r[i + 1] = (int)(i + 1) - pi[i];
    return r;
}

// ---- 回転・重なり ----
template <class A, class B>
std::vector<int> rotation_offsets(const A& a, const B& b) {  // a を左に k 回まわす（a[k..] + a[..k]）と b になる k（0 ≤ k < |a|。昇順。a, b がともに空なら {0}）  O(N)
    const auto x = str_internal::codes_of(a), y = str_internal::codes_of(b);
    std::vector<int> r;
    if (x.size() != y.size()) return r;
    if (x.empty()) return {0};
    std::vector<long long> xx(x);
    xx.insert(xx.end(), x.begin(), x.end() - 1);
    return str_internal::kmp_find(xx, y);
}
template <class A, class B>
bool is_rotation(const A& a, const B& b) { return !rotation_offsets(a, b).empty(); }  // 回転して一致するか
template <class T, class P>
std::vector<int> find_all_cyclic(const T& text, const P& pat) {  // text を無限に繰り返した列で pat が始まる位置 i（0 ≤ i < |text|。|pat| > |text| でもよい）  O(N + M)
    const auto t = str_internal::codes_of(text), p = str_internal::codes_of(pat);
    const std::size_t n = t.size(), m = p.size();
    if (n == 0) return {};
    std::vector<long long> tt;
    tt.reserve(n + m);
    while (tt.size() < n + (m > 0 ? m - 1 : 0)) tt.insert(tt.end(), t.begin(), t.end());
    std::vector<int> r;
    for (int i : str_internal::kmp_find(tt, p))
        if ((std::size_t)i < n) r.push_back(i);
    return r;
}
template <class A, class B>
int max_overlap(const A& a, const B& b) {  // a の接尾辞と b の接頭辞が一致する最大の長さ（a + b を短くつなげるとき）  O(N + M)
    const auto x = str_internal::codes_of(a), y = str_internal::codes_of(b);
    const int m = (int)y.size();
    if (m == 0) return 0;
    const std::vector<int> pi = str_internal::pi_codes(y);
    int k = 0;
    for (long long c : x) {
        if (k == m) k = pi[k - 1];
        while (k > 0 && c != y[k]) k = pi[k - 1];
        if (c == y[k]) ++k;
    }
    return k;
}

// ---- 置き換え・消すを、現れなくなるまで繰り返す（いちばん左の出現から） ----
//   スタックに KMP の状態を積み、置き換えた to は読み直す。O(N + 置き換えの回数 × (|from| + |to|))
//   |to| < |from| なら必ず終わる。置き換えが起きて to が from を含むと終わらないので止める（伸び続けるものは長さ 2×10^7 で止める）
//   |to| ≥ |from| で置き換えの回数が多いと遅い（ba → ab は Θ(N^2) 回）。to の値は text の要素の型に収まること
template <class T, class F, class U>
str_internal::seq_t<T> rewrite_until_stable(const T& text0, const F& from0, const U& to0) {  // from を to に（ABC394 C：WA → AC）
    using R = str_internal::seq_t<T>;
    using E = typename R::value_type;
    const auto& text = str_internal::seq(text0);
    const std::vector<long long> f = str_internal::codes_of(from0), to = str_internal::codes_of(to0);
    const int m = (int)f.size();
    if (m == 0) str_internal::fail("rewrite_until_stable", "from が空");
    const bool loops = !str_internal::kmp_find(to, f).empty();  // 1 回でも置き換えると終わらない
    const std::vector<int> pi = str_internal::pi_codes(f);
    std::vector<E> out;
    std::vector<int> st;          // st[k] = out[0, k] を読んだあとの KMP の状態
    std::vector<long long> back;  // 読み直す値（後ろから取り出す）
    out.reserve(text.size()), st.reserve(text.size());
    std::size_t idx = 0;
    while (idx < text.size() || !back.empty()) {
        E e;
        long long c;
        if (!back.empty()) c = back.back(), back.pop_back(), e = (E)c;
        else e = text[idx++], c = str_internal::code_of(e);
        int q = st.empty() ? 0 : st.back();
        while (q > 0 && f[q] != c) q = pi[q - 1];
        if (f[q] == c) ++q;
        out.push_back(e), st.push_back(q);
        if (q == m) {  // いちばん左の出現（スタックの中には出現がない）
            if (loops) str_internal::fail("rewrite_until_stable", "to が from を含む（いつまでも終わらない）");
            out.resize(out.size() - m), st.resize(st.size() - m);
            back.insert(back.end(), to.rbegin(), to.rend());
            if (out.size() + back.size() + (text.size() - idx) > 20000000) str_internal::fail("rewrite_until_stable", "長さが 2×10^7 を超えた（終わらない置き換え？）");
        }
    }
    return R(out.begin(), out.end());
}
template <class T, class P>
str_internal::seq_t<T> erase_until_stable(const T& text, const P& pat) {  // pat を消すことを現れなくなるまで（ABC328 D：ABC を消す）  O(N + M)
    return rewrite_until_stable(text, pat, std::vector<long long>());
}

// ---- 出現の一覧から、区間 [l, r) に収まる出現の数・i 以降で最初の出現を O(log N) で（ABC477 C） ----
struct Occurrences {
    Occurrences() = default;
    template <class T, class P>
    Occurrences(const T& text, const P& pat) : pos_(find_all(text, pat)), n_((int)str_internal::seq(text).size()), m_((int)str_internal::seq(pat).size()) {}
    const std::vector<int>& positions() const { return pos_; }  // 開始位置（昇順）
    int size() const { return (int)pos_.size(); }               // 出現の数
    long long count(int l, int r) const {                        // [l, r) に収まる（開始が l 以上・終わりが r 以下）出現の数  O(log N)
        str_internal::check_range("Occurrences::count", l, r, n_);
        if (r - l < m_) return 0;
        return idx(r - m_ + 1) - idx(l);
    }
    bool exists(int l, int r) const { return count(l, r) > 0; }  // [l, r) に収まる出現があるか
    int next(int i) const {                                       // i 以上で最初の出現の開始位置（なければ -1。0 ≤ i ≤ N）
        check_pos("Occurrences::next", i);
        const int k = idx(i);
        return k < size() ? pos_[k] : -1;
    }
    int prev(int i) const {  // i より前（開始位置 < i）で最後の出現の開始位置（なければ -1。0 ≤ i ≤ N）
        check_pos("Occurrences::prev", i);
        const int k = idx(i);
        return k > 0 ? pos_[k - 1] : -1;
    }

   private:
    int idx(int x) const { return str_internal::lower_index(pos_.data(), 0, size(), x); }  // 開始位置が x 以上の最初の番号
    void check_pos(const char* func, int i) const {
        if (i < 0 || i > n_) str_internal::fail(func, "位置 " + std::to_string(i) + " は 0〜N（= " + std::to_string(n_) + "）の外");
    }
    std::vector<int> pos_;
    int n_ = 0, m_ = 0;
};

// ---- KMP オートマトン：状態 = pat の先頭何文字まで一致しているか（0〜|pat|）。|pat| に着いたら pat が現れた ----
struct KmpAutomaton {
    KmpAutomaton() = default;
    template <class P>
    explicit KmpAutomaton(const P& pat0, const Alphabet& A = Alphabet::lower()) : A_(A) {
        const auto& pat = str_internal::seq(pat0);
        m_ = (int)pat.size();
        sg_ = A_.size();
        std::vector<int> p(m_);
        for (int i = 0; i < m_; ++i) p[i] = A_.id(pat[i]);
        nx_.assign((std::size_t)(m_ + 1) * sg_, 0);
        int fail = 0;  // 状態 q の「1 つ短い一致」
        for (int q = 0; q <= m_; ++q) {
            for (int c = 0; c < sg_; ++c) nx_[(std::size_t)q * sg_ + c] = (q < m_ && p[q] == c) ? q + 1 : (q == 0 ? 0 : nx_[(std::size_t)fail * sg_ + c]);
            if (q > 0 && q < m_) fail = nx_[(std::size_t)fail * sg_ + p[q]];
        }
    }
    int size() const { return m_ + 1; }  // 状態の数（0〜|pat|）
    template <class C>
    int next(int state, C c) const { return next_id(state, A_.id(c)); }  // 状態 state で文字 c を読んだあとの状態  O(1)
    int next_id(int state, int id) const {
        if (state < 0 || state > m_ || id < 0 || id >= sg_) str_internal::fail("KmpAutomaton::next", "状態か文字の番号が範囲外");
        return nx_[(std::size_t)state * sg_ + id];
    }
    bool accepting(int state) const { return state == m_; }  // pat がちょうど現れたところか
    const Alphabet& alphabet() const { return A_; }

   private:
    Alphabet A_;
    int m_ = 0, sg_ = 0;
    std::vector<int> nx_;
};
#endif  // GG_PATTERN_MATCH
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト
// =============================================================================
namespace pm_test {

#define PM_CHECK(...)                                                                          \
    do {                                                                                        \
        if (!(__VA_ARGS__)) {                                                                   \
            std::cerr << "[pattern_match] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
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

static std::vector<int> bf_find(const std::string& t, const std::string& p) {
    std::vector<int> r;
    for (int i = 0; i + (int)p.size() <= (int)t.size(); ++i)
        if (t.compare(i, p.size(), p) == 0) r.push_back(i);
    return r;
}
static int bf_lcp(const std::string& a, int i, const std::string& b, int j) {
    int k = 0;
    while (i + k < (int)a.size() && j + k < (int)b.size() && a[i + k] == b[j + k]) ++k;
    return k;
}

static void test_basic(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 8000 * TEST_SCALE; ++it) {
        const int sg = 1 + (int)(rng() % 3);
        const std::string t = rnd_str_shaped(rng, (int)(rng() % 25), sg), p = rnd_str_shaped(rng, (int)(rng() % 5), sg);
        const int n = (int)t.size(), m = (int)p.size();
        // z / π
        const auto z = z_array(t);
        for (int i = 0; i < n; ++i) PM_CHECK(z[i] == bf_lcp(t, 0, t, i));
        const auto pi = prefix_function(t);
        for (int i = 0; i < n; ++i) {
            int b = 0;
            for (int k = 1; k <= i; ++k)
                if (t.compare(0, k, t, i + 1 - k, k) == 0) b = k;
            PM_CHECK(pi[i] == b);
        }
        // 一致の長さ・出現
        const auto ml = match_lengths(t, p);
        for (int i = 0; i < n; ++i) PM_CHECK(ml[i] == bf_lcp(t, i, p, 0));
        std::vector<int> ex = bf_find(t, p);
        if (m == 0) {
            ex.clear();
            for (int i = 0; i <= n; ++i) ex.push_back(i);
        }
        PM_CHECK(find_all(t, p) == ex && count_occurrences(t, p) == (long long)ex.size());
        PM_CHECK(find_last(t, p) == (ex.empty() ? -1 : ex.back()));
        const int from = (int)(rng() % (n + 1));
        int ff = -1;
        for (int x : ex)
            if (x >= from) { ff = x; break; }
        PM_CHECK(find_first(t, p, from) == ff && find_first(t, p) == (ex.empty() ? -1 : ex[0]));
        if (m > 0) {
            long long c = 0;
            std::string rep;
            for (int i = 0; i < n;) {
                if (t.compare(i, m, p) == 0 && i + m <= n) ++c, rep += "XY", i += m;
                else rep += t[i++];
            }
            PM_CHECK(count_non_overlapping(t, p) == c && replace_all(t, p, "XY") == rep);
        }
        // 周期・境界
        std::vector<int> bd, pd;
        for (int k = n - 1; k >= 1; --k)
            if (t.compare(0, k, t, n - k, k) == 0) bd.push_back(k);
        for (int q = 1; q <= n; ++q) {
            bool ok = true;
            for (int i = 0; i + q < n; ++i) ok &= t[i] == t[i + q];
            if (ok) pd.push_back(q);
        }
        PM_CHECK(all_borders(t) == bd && all_periods(t) == pd && min_period(t) == (n ? pd[0] : 0));
        int root = n;
        for (int q : pd)
            if (n % q == 0) { root = q; break; }
        PM_CHECK(repetition_root(t) == root);
        const auto pmp = prefix_min_periods(t);
        PM_CHECK((int)pmp.size() == n + 1 && pmp[0] == 0);
        for (int i = 1; i <= n; i += 3) PM_CHECK(pmp[i] == min_period(t.substr(0, i)));
        // 回転
        std::string u = t;
        if (n && rng() % 2) std::rotate(u.begin(), u.begin() + rng() % n, u.end());
        else if (n && rng() % 2) u[rng() % n] = (char)('a' + rng() % 3);
        std::vector<int> ro;
        for (int k = 0; k < n; ++k)
            if (t.substr(k) + t.substr(0, k) == u) ro.push_back(k);
        if (n == 0) ro = {0};
        PM_CHECK(rotation_offsets(t, u) == ro && is_rotation(t, u) == !ro.empty());
        // 無限に繰り返した中の出現（pat は長くてもよい）
        const std::string q = rnd_str_shaped(rng, (int)(rng() % 12), sg);
        if (n > 0) {
            std::string big;
            while (big.size() < (size_t)n + q.size()) big += t;
            std::vector<int> cy;
            for (int i = 0; i < n; ++i)
                if (big.compare(i, q.size(), q) == 0) cy.push_back(i);
            PM_CHECK(find_all_cyclic(t, q) == cy);
        }
        // 重なり
        int ov = 0;
        for (int k = 1; k <= std::min(n, (int)q.size()); ++k)
            if (t.compare(n - k, k, q, 0, k) == 0) ov = k;
        PM_CHECK(max_overlap(t, q) == ov);
        ++cases;
    }
}

static void test_automaton(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 3000 * TEST_SCALE; ++it) {
        const int sg = 1 + (int)(rng() % 3);
        const std::string p = rnd_str_shaped(rng, (int)(rng() % 6), sg), t = rnd_str(rng, (int)(rng() % 30), sg);
        const KmpAutomaton KA(p);
        PM_CHECK(KA.size() == (int)p.size() + 1);
        int st = 0;
        long long acc = 0;
        for (int i = 0; i < (int)t.size(); ++i) {
            st = KA.next(st, t[i]);
            acc += KA.accepting(st);
            // 状態 = 読んだ文字列の接尾辞で p の接頭辞であるものの最長（p の全体まで）
            int best = 0;
            for (int k = 1; k <= std::min((int)p.size(), i + 1); ++k)
                if (t.compare(i + 1 - k, k, p, 0, k) == 0) best = k;
            PM_CHECK(st == best);
        }
        if (!p.empty()) PM_CHECK(acc == count_occurrences(t, p));
        ++cases;
    }
    // 括弧の文字で
    const KmpAutomaton KB("()", Alphabet::of("()"));
    PM_CHECK(KB.next(KB.next(0, '('), ')') == 2 && KB.accepting(2));
}

static void test_values(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 3000 * TEST_SCALE; ++it) {
        const int n = (int)(rng() % 20), m = 1 + (int)(rng() % 3);
        std::vector<long long> t(n), p(m);
        const long long big = 1000000000000000000LL;
        for (auto& x : t) x = (rng() % 2 ? big : -big) * (long long)(rng() % 2);
        for (auto& x : p) x = (rng() % 2 ? big : -big) * (long long)(rng() % 2);
        std::vector<int> ex;
        for (int i = 0; i + m <= n; ++i)
            if (std::equal(p.begin(), p.end(), t.begin() + i)) ex.push_back(i);
        PM_CHECK(find_all(t, p) == ex);
        ++cases;
    }
    // string と vector<int> を混ぜる・リテラル
    PM_CHECK(find_all(std::vector<int>{'a', 'b', 'a', 'b'}, "ab") == std::vector<int>({0, 2}));
    PM_CHECK(count_occurrences("aaaa", "aa") == 3 && count_non_overlapping("aaaa", "aa") == 2 && replace_all("aaaa", "aa", "b") == "bb");
    PM_CHECK(replace_all(std::vector<int>{1, 2, 1, 2}, std::vector<int>{1, 2}, std::vector<int>{9}) == std::vector<int>({9, 9}));
}

// 最終照合で足したもの：置き換えの繰り返し・区間の出現
static std::string bf_rewrite(std::string s, const std::string& f, const std::string& t) {
    while (true) {
        const std::size_t i = s.find(f);
        if (i == std::string::npos) return s;
        s.replace(i, f.size(), t);
    }
}
static void test_final_additions(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 4000 * TEST_SCALE; ++it) {
        const int sg = 1 + (int)(rng() % 3);
        const std::string s = rnd_str_shaped(rng, (int)(rng() % 30), sg);
        const std::string f = rnd_str_shaped(rng, 1 + (int)(rng() % 3), sg);
        std::string t = rnd_str(rng, (int)(rng() % (f.size() + 1)), sg);  // |t| ≤ |f|
        if (t.find(f) != std::string::npos) t = "";
        if (t.size() == f.size() && rng() % 2) t = std::string(f.rbegin(), f.rend());  // 同じ長さ（ab → ba など）
        if (t.find(f) != std::string::npos) t = "";
        PM_CHECK(rewrite_until_stable(s, f, t) == bf_rewrite(s, f, t));
        PM_CHECK(erase_until_stable(s, f) == bf_rewrite(s, f, ""));
        // 区間の出現
        const std::string p = rnd_str_shaped(rng, (int)(rng() % 3), sg);
        const Occurrences occ(s, p);
        const int n = (int)s.size();
        PM_CHECK(occ.positions() == find_all(s, p) && occ.size() == (int)find_all(s, p).size());
        for (int q = 0; q < 20; ++q) {
            int l = (int)(rng() % (n + 1)), r = (int)(rng() % (n + 1));
            if (l > r) std::swap(l, r);
            long long c = 0;
            for (int i = l; i + (int)p.size() <= r; ++i) c += s.compare(i, p.size(), p) == 0;
            PM_CHECK(occ.count(l, r) == c && occ.exists(l, r) == (c > 0));
            int nx = -1, pv = -1;
            for (int x : find_all(s, p)) {
                if (x >= l && nx < 0) nx = x;
                if (x < l) pv = x;
            }
            PM_CHECK(occ.next(l) == nx && occ.prev(l) == pv);
        }
        ++cases;
    }
    // 数の列・1 文字の pat・伸びる置き換え
    PM_CHECK(rewrite_until_stable(std::vector<int>{1, 2, 1, 2, 3}, std::vector<int>{1, 2}, std::vector<int>{2}) == std::vector<int>({2, 2, 3}));
    PM_CHECK(erase_until_stable(std::string("AABCBCC"), "ABC") == "C" && rewrite_until_stable(std::string("WWWWA"), "WA", "AC") == "ACCCC");
    PM_CHECK(find_all("abca", 'a') == std::vector<int>({0, 3}) && replace_all("a,b", ',', '.') == "a.b");
    PM_CHECK(rewrite_until_stable(std::string("xyz"), "a", "aa") == "xyz");  // 置き換えが起きなければ to が from を含んでもよい
    ++cases;
}

static void random_test() {
    std::mt19937_64 rng(20261012);
    int cases = 0;
    test_basic(rng, cases);
    test_automaton(rng, cases);
    test_values(rng, cases);
    test_final_additions(rng, cases);
    std::cout << "[pattern_match] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[pattern_match] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    auto ms = [](std::chrono::steady_clock::time_point a) {
        return (long long)std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count();
    };
    const int n = 10000000;
    const std::string s = rnd_str_shaped(rng, n, 2), a(n, 'a');
    auto t0 = std::chrono::steady_clock::now();
    const auto z = z_array(s);
    const auto z2 = z_array(a);
    const long long t_z = ms(t0);
    t0 = std::chrono::steady_clock::now();
    const auto pi = prefix_function(s);
    const long long t_pi = ms(t0);
    t0 = std::chrono::steady_clock::now();
    const auto occ = find_all(s, s.substr(1234567, 20));
    const long long t_find = ms(t0);
    t0 = std::chrono::steady_clock::now();
    const auto ml = match_lengths(s, s.substr(0, 100000));
    const long long t_ml = ms(t0);
    const std::string p = rnd_str(rng, 100000, 26);
    t0 = std::chrono::steady_clock::now();
    const KmpAutomaton KA(p);
    int st = 0;
    for (int i = 0; i < n; ++i) st = KA.next_id(st, (int)(rng() % 26));
    const long long t_ka = ms(t0);
    std::string nest = std::string(500000, '(') + "xx" + std::string(500000, ')'), wa = std::string(1000000, 'W') + "A";
    t0 = std::chrono::steady_clock::now();
    const auto r1 = rewrite_until_stable(nest, "(xx)", "xx"), r2 = rewrite_until_stable(wa, "WA", "AC");
    const long long t_rw = ms(t0);
    std::cout << "[pattern_match] benchmark: z_array（N = 10^7 × 2）" << t_z << " ms, prefix_function（N = 10^7）" << t_pi << " ms, find_all（N = 10^7）"
              << t_find << " ms, match_lengths（N = 10^7）" << t_ml << " ms, KmpAutomaton（M = 10^5, σ = 26）+ 10^7 回の遷移 " << t_ka
              << " ms, rewrite_until_stable（入れ子 10^6・WA→AC 10^6）" << t_rw << " ms  (check "
              << (z[5] + z2[3] + pi[7] + (long long)occ.size() + ml[9] + st + (long long)r1.size() + (long long)r2.size()) % 7 << ")\n";
}

}  // namespace pm_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace pm_demo {

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
#define SHOW(expr, comment) pm_demo::show(#expr, (expr), comment)
#define RUN(...)                   \
    pm_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("パターン検索（pattern_match） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  z_array / prefix_function / match_lengths              Z 配列 / KMP の境界 / 各位置で pat と何文字一致するか\n"
              << "  find_all / count_occurrences / find_first / find_last  出現位置・回数（重なりも）\n"
              << "  count_non_overlapping / replace_all                    重ならない出現の数 / 全部置き換え（1 回）\n"
              << "  rewrite_until_stable / erase_until_stable              置き換え・消すを現れなくなるまで繰り返す\n"
              << "  Occurrences                                            区間 [l, r) に収まる出現の数・次の出現\n"
              << "  all_borders / min_period / all_periods / repetition_root / prefix_min_periods   境界・周期\n"
              << "  rotation_offsets / is_rotation / find_all_cyclic        回転・無限に繰り返した中の出現\n"
              << "  max_overlap                                             a の後ろと b の前の重なり\n"
              << "  KmpAutomaton                                            パターンを含まない文字列の DP\n";

    api("【1】 z_array(s) / prefix_function(s)", "vector<int>",
        "z[i] = s と s[i..) の最長共通接頭辞（z[0] = |s|。ACL の z_algorithm と同じ）/ π[i] = s[0, i] の最長の真の境界", "O(N)");
    SHOW(z_array("aabaab"), "");
    SHOW(prefix_function("aabaab"), "");

    api("【2】 match_lengths(text, pat)", "vector<int>（長さ |text|）",
        "各 i で text[i..) と pat の最長共通接頭辞の長さ（Z を text と pat に）", "O(N + M)");
    SHOW(match_lengths("abcab", "ab"), "");
    {
        // ABC257 G: T を S の接頭辞をつなげて作る最小の個数（位置 i から i + 1〜i + m[i] へ進める → 貪欲）
        const std::string S = "aba", T = "ababaab";
        const auto m = match_lengths(T, S);
        int cnt = 0, cur_end = 0, far = 0;
        bool ok = true;
        for (int i = 0; i < (int)T.size() && ok; ++i) {
            far = std::max(far, i + m[i]);
            if (i == cur_end) {
                if (far == i) ok = false;
                ++cnt, cur_end = far;
            }
        }
        note("ABC257 G（Prefix Concatenation）: m = match_lengths(T, S) で「i から i + m[i] まで進める」→ 貪欲に遠くへ");
        std::cout << "    S = aba, T = ababaab → " << (ok ? cnt : -1) << "\n";
    }

    api("【3】 find_all / count_occurrences / find_first / find_last / count_non_overlapping / replace_all", "vector<int> / long long / int / string",
        "pat の出現（重なりも）。重ならない数・全部置き換えは左から貪欲に。数の列にも（ABC150 F のように差や XOR の列にして）", "O(N + M)");
    SHOW(find_all("abababa", "aba"), "重なりも");
    SHOW(count_occurrences("abababa", "aba"), "");
    SHOW(count_non_overlapping("abababa", "aba"), "左から重ならないように");
    SHOW(find_first("abababa", "aba", 1), "1 以上で最初");
    SHOW(find_last("abababa", "aba"), "");
    SHOW(find_first("abc", "x"), "ない → -1");
    SHOW(replace_all("abababa", "aba", "X"), "");
    SHOW(find_all(std::vector<int>{3, 1, 4, 1, 4}, std::vector<int>{1, 4}), "数の列");
    note("ABC343 G のように「a が b に含まれるか」は find_first(b, a) >= 0。1 文字の pat も渡せる（find_all(s, 'a')）。");

    api("【3b】 rewrite_until_stable(text, from, to) / erase_until_stable(text, pat)", "string（vector なら vector）",
        "いちばん左の from を to にすることを、from が現れなくなるまで繰り返した結果（replace_all をループで回すと O(N^2) になる形も速い）",
        "O(N + 置き換えの回数 × (|from| + |to|))");
    SHOW(rewrite_until_stable("WWWWA", "WA", "AC"), "ABC394 C（Debug）：WA → AC をなくなるまで");
    SHOW(erase_until_stable("BAABCBCCABCAC", "ABC"), "ABC328 D（Take ABC）の入力例 1：ABC を消し続ける");
    SHOW(rewrite_until_stable("((xx))(xx)", "(xx)", "xx"), "ABC454 D：(xx) → xx を繰り返した形どうしを比べる");
    note("|to| < |from| なら必ず終わる。to が from を含むと終わらないので止まる。");

    api("【3c】 Occurrences occ(text, pat)：count(l, r) / exists(l, r) / next(i) / prev(i) / positions()", "構造体",
        "区間 [l, r) に収まる（開始が l 以上・終わりが r 以下）出現の数・あるか / i 以上で最初・i より前で最後の出現の開始位置", "前処理 O(N + M)・O(log N)");
    RUN(const Occurrences occ("abcabcab", "ab"));
    SHOW(occ.count(1, 8), "[1, 8) に収まる \"ab\"（3 と 6）");
    SHOW(occ.exists(1, 4), "\"bca\" には \"ab\" がない");
    SHOW(occ.next(1), "1 以上で最初の出現");
    note("ABC477 C（区間に T が現れるか）。部分文字列 s[l, r) が全体に何回現れるかは 03_suffix_array の count_occurrences(l, r)（意味が違う）。");

    api("【4】 all_borders / min_period / all_periods / repetition_root / prefix_min_periods", "vector<int> / int",
        "境界（接頭辞 = 接尾辞）・周期（s[i] = s[i + p]）・s = t^k となる最短の t・各接頭辞の最小の周期", "O(N)");
    SHOW(all_borders("abaabaab"), "大きい順");
    SHOW(min_period("abaabaab"), "割り切れなくてよい");
    SHOW(all_periods("abaabaab"), "");
    SHOW(repetition_root("abcabcabc"), "\"abc\" × 3");
    SHOW(repetition_root("abcabca"), "繰り返しでない → 長さそのもの");
    SHOW(prefix_min_periods("abab"), "r[i] = s[0, i) の最小の周期");

    api("【5】 rotation_offsets(a, b) / is_rotation(a, b) / find_all_cyclic(text, pat)", "vector<int> / bool / vector<int>",
        "a を左に k 回まわす（a[k..] + a[..k]）と b になる k / text を無限に繰り返した中で pat が始まる位置（|pat| > |text| でもよい）", "O(N + M)");
    SHOW(rotation_offsets("abcab", "cabab"), "");
    SHOW(rotation_offsets("abab", "baba"), "周期があると複数");
    SHOW(is_rotation("abc", "acb"), "");
    SHOW(find_all_cyclic("ab", "babab"), "\"ababab…\" の 1 文字目から");
    note("ABC135 F（Strings of Eternity）: find_all_cyclic(s, t) で i → i + |t| (mod |s|) の辺を作り、");
    note("graph_grid_v2 の 06_functional_graph で閉路と最長の道を調べる。");

    api("【6】 max_overlap(a, b)", "int", "a の接尾辞と b の接頭辞が一致する最大の長さ（a と b を重ねて短くつなぐ）", "O(N + M)");
    SHOW(max_overlap("abcab", "cabd"), "\"cab\" が重なる");
    SHOW(max_overlap("abc", "xyz"), "");
    note("ABC343 G（Compress Strings）: 含まれる文字列を除き、max_overlap を重みにして bit DP。");

    api("【7】 KmpAutomaton KA(pat, A) / KA.next(state, c) / KA.accepting(state)", "構造体",
        "状態 = 読んだ文字列の末尾が pat の先頭何文字と一致しているか。pat を含まない文字列の数え上げ DP に", "構築 O(Mσ)・遷移 O(1)");
    {
        // 長さ 6 の a, b からなる文字列で "aba" を含まないものの数（DP）
        const KmpAutomaton KA("aba", Alphabet::of("ab"));
        std::vector<long long> dp(KA.size(), 0);
        dp[0] = 1;
        for (int step = 0; step < 6; ++step) {
            std::vector<long long> nd(KA.size(), 0);
            for (int st = 0; st < KA.size(); ++st)
                for (char c : {'a', 'b'}) {
                    const int to = KA.next(st, c);
                    if (!KA.accepting(to)) nd[to] += dp[st];
                }
            dp = nd;
        }
        note("長さ 6 の {a, b} の文字列で \"aba\" を含まないもの: dp[状態] を 6 回遷移（accepting に入る遷移を捨てる）");
        std::cout << "    → " << std::accumulate(dp.begin(), dp.end(), 0LL) << "（全 64 個から）\n";
        note("長さが 10^18 なら、遷移を行列にして math_v2 の 05_matrix（Matrix<mint>::pow）。パターンが多いなら 02_aho_corasick。");
    }
}

}  // namespace pm_demo

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
    pm_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    pm_test::benchmark();    // 速度計測
    pm_demo::demo();         // 使い方デモ
    return 0;
}
