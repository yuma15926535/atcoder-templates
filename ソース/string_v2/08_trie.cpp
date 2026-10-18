// =============================================================================
//  Trie（文字列の多重集合：追加・削除・接頭辞の数・辞書順 k 番目・順位・LCP の最大と和・ノードをたどる）   --- 文字列第二世代
// -----------------------------------------------------------------------------
//  分類       : 08 文字列の集合（接頭辞）
//  依存       : 00_str_base（同梱済み）
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）
//  文字コード : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  Trie TR(A)：A は文字の種類（Alphabet。既定は英小文字。辞書順は A の番号の順）。範囲外の文字はメッセージを出して止める
//  子は 1 次元の配列（ノード数 × σ 個の int。全部の文字列の長さの合計 10^6・σ = 26 で約 100MB）。ノード 0 が根（空文字列）
//  同じ文字列を何個でも入れられる（数を持つ）。消してもノードは残る（数が 0 になるだけ）
//
//  機能一覧                                                                       計算量
//    TR.insert(s, cnt) / erase(s, cnt)     cnt 個入れる（ノードを返す）/ 消す（足りなければ何もせず false）   O(|s|)
//    TR.erase_prefix(p)                    p で始まる文字列を全部消す（消した数）                          O(|p| + 消えたノード数 × σ)
//    TR.count(s) / contains(s)             s の数 / 入っているか                                         O(|s|)
//    TR.count_prefix(p)                    p で始まる文字列の数                                          O(|p|)
//    TR.count_prefixes_of(s)               s の接頭辞である文字列の数（s 自身も）                          O(|s|)
//    TR.size() / empty()                   入っている数（同じものも数える）
//    TR.kth(k) / kth_ids(k)                辞書順 k 番目（0 から）                                       O(|答え| σ)
//    TR.rank(s)                            s より小さい文字列の数                                         O(|s| σ)
//    TR.words() / words_ids()              全部を辞書順に                                                O(ノード数 × σ + 出力)
//    TR.max_lcp(s, exclude_self)           入っている文字列との LCP の最大（s を 1 つ除いても）             O(|s|)
//    TR.sum_lcp_with(s) / sum_lcp_pairs()  全部との LCP の和 / 全部の組の LCP の和                        O(|s|) / O(ノード数)
//    ノード: find(s) / path(s) / child(v, c) / child_by_id(v, id) / parent(v) / depth(v) / label(v)
//            pass_count(v) / end_count(v) / string_of(v) / preorder() / node_count() / root()
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    trie_test::random_test()  … map で持った多重集合と全機能を突き合わせるランダムテスト
//    trie_test::benchmark()    … 速度計測
//    trie_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
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
#ifndef GG_TRIE
#define GG_TRIE
// ---- Trie（文字列の多重集合。文字列第二世代） ----
// Trie TR(A)：A は文字の種類（既定は英小文字）。子は 1 次元の配列（ノード数 × σ 個の int。10^6 文字・σ = 26 で約 100MB）
// 同じ文字列を何個でも入れられる（数を持つ）。消してもノードは残る（数が 0 になるだけ）。ノード 0 が根（空文字列）
struct Trie {
    explicit Trie(const Alphabet& A = Alphabet::lower()) : A_(A), sg_(A.size()) { new_node(-1, -1); }

    // ---- 追加・削除 ----
    template <class Seq>
    int insert(const Seq& s0, long long cnt = 1) {  // s を cnt 個入れる（0 個でもノードは作る）。s のノードを返す  O(|s|)
        if (cnt < 0) str_internal::fail("Trie::insert", "cnt ≥ 0 で使う（減らすときは erase）");
        const auto& s = str_internal::seq(s0);
        int v = 0;
        nd_[0].pass += cnt;
        for (std::size_t i = 0; i < s.size(); ++i) {
            const int c = A_.id(s[i]);
            int u = next_[(std::size_t)v * sg_ + c];
            if (u < 0) u = new_node(v, c), next_[(std::size_t)v * sg_ + c] = u;
            v = u;
            nd_[v].pass += cnt;
        }
        nd_[v].end += cnt;
        return v;
    }
    template <class Seq>
    bool erase(const Seq& s, long long cnt = 1) {  // s を cnt 個消す（入っている数が cnt より少なければ何もせず false）  O(|s|)
        if (cnt < 0) str_internal::fail("Trie::erase", "cnt ≥ 0 で使う");
        const int v = find(s);
        if (v < 0) return cnt == 0;
        if (nd_[v].end < cnt) return false;
        nd_[v].end -= cnt;
        for (int u = v; u >= 0; u = nd_[u].par) nd_[u].pass -= cnt;
        return true;
    }
    template <class Seq>
    long long erase_prefix(const Seq& p) {  // p で始まる文字列を全部消す。消した数を返す  O(|p| + 消えたノード数 × σ)
        const int v = find(p);
        if (v < 0 || nd_[v].pass == 0) return 0;
        const long long c = nd_[v].pass;
        for (int u = nd_[v].par; u >= 0; u = nd_[u].par) nd_[u].pass -= c;
        std::vector<int> st = {v};  // 数が残っている部分木だけを 0 にする
        while (!st.empty()) {
            const int u = st.back();
            st.pop_back();
            nd_[u].pass = nd_[u].end = 0;
            for (int k = 0; k < sg_; ++k) {
                const int w = next_[(std::size_t)u * sg_ + k];
                if (w >= 0 && nd_[w].pass > 0) st.push_back(w);
            }
        }
        return c;
    }
    void clear() { *this = Trie(A_); }

    // ---- 数える ----
    template <class Seq>
    long long count(const Seq& s) const {  // s の数  O(|s|)
        const int v = find(s);
        return v < 0 ? 0 : nd_[v].end;
    }
    template <class Seq>
    bool contains(const Seq& s) const { return count(s) > 0; }
    template <class Seq>
    long long count_prefix(const Seq& p) const {  // p で始まる文字列の数（p 自身も）  O(|p|)
        const int v = find(p);
        return v < 0 ? 0 : nd_[v].pass;
    }
    template <class Seq>
    long long count_prefixes_of(const Seq& s) const {  // s の接頭辞である文字列の数（空文字列・s 自身も）  O(|s|)
        long long r = 0;
        for (int v : path(s)) r += nd_[v].end;
        return r;
    }
    long long size() const { return nd_[0].pass; }  // 入っている文字列の数（同じものも数える）
    bool empty() const { return nd_[0].pass == 0; }

    // ---- 辞書順 ----
    std::vector<int> kth_ids(long long k) const {  // 辞書順 k 番目（0 から。同じものも数える）を番号の列で  O(|答え| σ)
        if (k < 0 || k >= size()) str_internal::fail("Trie::kth", "k = " + std::to_string(k) + " は 0〜size()-1 の外");
        std::vector<int> r;
        int v = 0;
        while (true) {
            if (k < nd_[v].end) return r;
            k -= nd_[v].end;
            for (int c = 0; c < sg_; ++c) {
                const int u = next_[(std::size_t)v * sg_ + c];
                if (u < 0) continue;
                if (k < nd_[u].pass) {
                    v = u, r.push_back(c);
                    break;
                }
                k -= nd_[u].pass;
            }
        }
    }
    std::string kth(long long k) const { return from_ids(kth_ids(k), A_); }  // 辞書順 k 番目の文字列
    template <class Seq>
    long long rank(const Seq& s0) const {  // s より辞書順で小さい文字列の数  O(|s| σ)
        const auto& s = str_internal::seq(s0);
        long long r = 0;
        int v = 0;
        for (std::size_t i = 0; i < s.size(); ++i) {
            r += nd_[v].end;  // s の真の接頭辞は s より小さい
            const int c = A_.id(s[i]);
            for (int k = 0; k < c; ++k) {
                const int u = next_[(std::size_t)v * sg_ + k];
                if (u >= 0) r += nd_[u].pass;
            }
            v = next_[(std::size_t)v * sg_ + c];
            if (v < 0) break;
        }
        return r;
    }
    std::vector<std::vector<int>> words_ids() const {  // 全部を辞書順に（番号の列。同じものは数だけ）
        std::vector<std::vector<int>> r;
        std::vector<int> cur;
        for (int v : preorder()) {
            if (nd_[v].pass == 0) continue;
            cur.resize(nd_[v].dep);
            if (v != 0) cur[nd_[v].dep - 1] = nd_[v].lab;
            for (long long t = 0; t < nd_[v].end; ++t) r.push_back(cur);
        }
        return r;
    }
    std::vector<std::string> words() const {  // 全部を辞書順に（文字列。同じものは数だけ）
        std::vector<std::string> r;
        for (const auto& w : words_ids()) r.push_back(from_ids(w, A_));
        return r;
    }

    // ---- 共通接頭辞（LCP） ----
    template <class Seq>
    int max_lcp(const Seq& s0, bool exclude_self = false) const {  // 入っている文字列との LCP の最大（exclude_self で s を 1 つ除いて。ABC287 E）  O(|s|)
        const auto& s = str_internal::seq(s0);
        if (exclude_self && count(s) == 0) str_internal::fail("Trie::max_lcp", "exclude_self = true は s が入っているときだけ");
        const long long sub = exclude_self ? 1 : 0;  // s の道の上のノードは数を 1 減らして見る
        int v = 0, d = 0;
        for (std::size_t i = 0; i < s.size(); ++i) {
            const int u = next_[(std::size_t)v * sg_ + A_.id(s[i])];
            if (u < 0 || nd_[u].pass - sub <= 0) break;
            v = u, ++d;
        }
        return d;
    }
    template <class Seq>
    long long sum_lcp_with(const Seq& s) const {  // 入っている全部の文字列 t との LCP(s, t) の和  O(|s|)
        long long r = 0;
        for (int v : path(s))
            if (v != 0) r = str_internal::checked_add(r, nd_[v].pass, "Trie::sum_lcp_with");
        return r;
    }
    long long sum_lcp_pairs() const {  // 入っている文字列の全部の組（i < j）の LCP の和（ABC353 E）  O(ノード数)
        long long r = 0;
        for (int v = 1; v < node_count(); ++v) {
            long long x = nd_[v].pass, y = x - 1, t;  // C(pass, 2)。先に偶数の方を 2 で割る
            if (x % 2 == 0) x /= 2;
            else y /= 2;
            if (__builtin_mul_overflow(x, y, &t)) str_internal::fail("Trie::sum_lcp_pairs", "値が long long に収まらない");
            r = str_internal::checked_add(r, t, "Trie::sum_lcp_pairs");
        }
        return r;
    }

    // ---- ノード ----
    int root() const { return 0; }
    int node_count() const { return (int)nd_.size(); }
    template <class Seq>
    int find(const Seq& s0) const {  // s のノード（ノードがなければ -1。数が 0 でもノードがあれば返す）  O(|s|)
        const auto& s = str_internal::seq(s0);
        int v = 0;
        for (std::size_t i = 0; i < s.size() && v >= 0; ++i) v = next_[(std::size_t)v * sg_ + A_.id(s[i])];
        return v;
    }
    template <class Seq>
    std::vector<int> path(const Seq& s0) const {  // 根から s をたどったノードの列（根を含む。途中でなくなればそこまで）  O(|s|)
        const auto& s = str_internal::seq(s0);
        std::vector<int> r = {0};
        int v = 0;
        for (std::size_t i = 0; i < s.size(); ++i) {
            v = next_[(std::size_t)v * sg_ + A_.id(s[i])];
            if (v < 0) break;
            r.push_back(v);
        }
        return r;
    }
    template <class C>
    int child(int v, C c) const {  // 子（文字 c。なければ -1）  O(1)
        check_node(v, "Trie::child");
        return next_[(std::size_t)v * sg_ + A_.id(c)];
    }
    int child_by_id(int v, int id) const {  // 子（文字の番号 id。なければ -1）  O(1)
        check_node(v, "Trie::child_by_id");
        if (id < 0 || id >= sg_) str_internal::fail("Trie::child_by_id", "番号は 0〜σ-1");
        return next_[(std::size_t)v * sg_ + id];
    }
    int parent(int v) const { return check_node(v, "Trie::parent"), nd_[v].par; }  // 親（根は -1）
    int depth(int v) const { return check_node(v, "Trie::depth"), nd_[v].dep; }    // 深さ = 文字列の長さ
    int label(int v) const { return check_node(v, "Trie::label"), nd_[v].lab; }    // 親からの文字の番号（根は -1）
    long long pass_count(int v) const { return check_node(v, "Trie::pass_count"), nd_[v].pass; }  // v を通る（v の文字列で始まる）文字列の数
    long long end_count(int v) const { return check_node(v, "Trie::end_count"), nd_[v].end; }     // v で終わる文字列の数
    std::string string_of(int v) const {  // v の文字列  O(深さ)
        check_node(v, "Trie::string_of");
        std::vector<int> ids;
        for (int u = v; u > 0; u = nd_[u].par) ids.push_back(nd_[u].lab);
        std::reverse(ids.begin(), ids.end());
        return from_ids(ids, A_);
    }
    std::vector<int> preorder() const {  // 全部のノードを辞書順の行きがけ順に（親が子より先。逆順にすれば子が先。数が 0 のノードも）  O(ノード数 × σ)
        std::vector<int> r, st = {0};
        r.reserve(nd_.size());
        while (!st.empty()) {
            const int v = st.back();
            st.pop_back();
            r.push_back(v);
            for (int c = sg_ - 1; c >= 0; --c) {
                const int u = next_[(std::size_t)v * sg_ + c];
                if (u >= 0) st.push_back(u);
            }
        }
        return r;
    }
    const Alphabet& alphabet() const { return A_; }

   private:
    int new_node(int p, int c) {
        const int v = (int)nd_.size();
        next_.insert(next_.end(), (std::size_t)sg_, -1);
        nd_.push_back({0, 0, p, p < 0 ? 0 : nd_[p].dep + 1, c});
        return v;
    }
    void check_node(int v, const char* func) const {
        if (v < 0 || v >= node_count()) str_internal::fail(func, "ノード " + std::to_string(v) + " は 0〜node_count()-1 の外");
    }
    Alphabet A_;
    int sg_;
    struct Node {
        long long pass, end;  // 通る数・終わる数
        int par, dep, lab;    // 親・深さ・親からの文字の番号
    };
    std::vector<int> next_;  // 子（ノード数 × σ）
    std::vector<Node> nd_;
};
#endif  // GG_TRIE
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト
// =============================================================================
namespace trie_test {

#define TRIE_CHECK(...)                                                                          \
    do {                                                                                        \
        if (!(__VA_ARGS__)) {                                                                   \
            std::cerr << "[trie] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
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

using Ids = std::vector<int>;
static int lcp_ids(const Ids& a, const Ids& b) {
    std::size_t k = 0;
    while (k < a.size() && k < b.size() && a[k] == b[k]) ++k;
    return (int)k;
}
static bool is_prefix(const Ids& p, const Ids& s) { return p.size() <= s.size() && std::equal(p.begin(), p.end(), s.begin()); }

// 愚直：番号の列 → 数 の map（map の順 = 番号の辞書順）
template <class Make>
static void run_ops(std::mt19937_64& rng, const Alphabet& A, int sg, Make make, int& cases) {
    Trie TR(A);
    std::map<Ids, long long> M;
    std::vector<Ids> pool;  // 使った文字列（何度も使う）
    auto rnd_ids = [&]() {
        if (!pool.empty() && rng() % 2) {
            Ids s = pool[rng() % pool.size()];
            if (rng() % 3 == 0 && !s.empty()) s.resize(rng() % s.size());  // 接頭辞
            else if (rng() % 3 == 0) s.push_back((int)(rng() % sg));
            return s;
        }
        Ids s(rng() % 6);
        for (auto& c : s) c = (int)(rng() % sg);
        pool.push_back(s);
        return s;
    };
    for (int q = 0; q < 300; ++q) {
        const int op = (int)(rng() % 10);
        const Ids s = rnd_ids();
        const auto S = make(s);  // 文字列 / 数の列
        if (op <= 2) {
            const long long cnt = rng() % 4 == 0 ? 0 : 1 + (long long)(rng() % 3);
            const int v = TR.insert(S, cnt);
            M[s] += cnt;
            TRIE_CHECK(v == TR.find(S) && TR.depth(v) == (int)s.size());
        } else if (op <= 4) {
            const long long cnt = (long long)(rng() % 3);
            const long long have = M.count(s) ? M[s] : 0;
            TRIE_CHECK(TR.erase(S, cnt) == (have >= cnt));
            if (have >= cnt) M[s] -= cnt;
        } else if (op == 5 && rng() % 3 == 0) {
            long long removed = 0;
            for (auto& [t, c] : M)
                if (is_prefix(s, t)) removed += c, c = 0;
            TRIE_CHECK(TR.erase_prefix(S) == removed);
        }
        // 数える
        long long total = 0, cnt_s = 0, pre = 0, pre_of = 0, rk = 0, sum_with = 0;
        int best = 0, best_ex = 0;
        bool have_s = false;
        for (const auto& [t, c] : M) {
            if (c == 0) continue;
            total += c;
            if (t == s) cnt_s = c, have_s = true;
            if (is_prefix(s, t)) pre += c;
            if (is_prefix(t, s)) pre_of += c;
            if (t < s) rk += c;
            sum_with += c * lcp_ids(s, t);
            best = std::max(best, lcp_ids(s, t));
            if (t != s || c >= 2) best_ex = std::max(best_ex, lcp_ids(s, t));
        }
        TRIE_CHECK(TR.size() == total && TR.empty() == (total == 0));
        TRIE_CHECK(TR.count(S) == cnt_s && TR.contains(S) == have_s);
        TRIE_CHECK(TR.count_prefix(S) == pre && TR.count_prefixes_of(S) == pre_of);
        TRIE_CHECK(TR.rank(S) == rk && TR.sum_lcp_with(S) == sum_with);
        TRIE_CHECK(TR.max_lcp(S) == best);
        if (have_s) TRIE_CHECK(TR.max_lcp(S, true) == best_ex);
        if (q % 10 == 0) {
            // 辞書順・組の LCP の和
            std::vector<Ids> all;
            for (const auto& [t, c] : M)
                for (long long k = 0; k < c; ++k) all.push_back(t);
            TRIE_CHECK(TR.words_ids() == all);
            for (long long k = 0; k < (long long)all.size(); ++k) TRIE_CHECK(TR.kth_ids(k) == all[k]);
            long long pairs = 0;
            for (std::size_t i = 0; i < all.size(); ++i)
                for (std::size_t j = i + 1; j < all.size(); ++j) pairs += lcp_ids(all[i], all[j]);
            TRIE_CHECK(TR.sum_lcp_pairs() == pairs);
            // ノード：行きがけ順は辞書順で、親が先
            const auto pre_order = TR.preorder();
            TRIE_CHECK((int)pre_order.size() == TR.node_count() && pre_order[0] == TR.root());
            std::vector<int> seen(TR.node_count(), 0);
            Ids prev_ids;
            for (std::size_t k = 0; k < pre_order.size(); ++k) {
                const int v = pre_order[k];
                TRIE_CHECK(!seen[v] && (v == 0 || seen[TR.parent(v)]));
                seen[v] = 1;
                Ids ids;
                for (int u = v; u > 0; u = TR.parent(u)) ids.push_back(TR.label(u));
                std::reverse(ids.begin(), ids.end());
                TRIE_CHECK((int)ids.size() == TR.depth(v) && (k == 0 || prev_ids < ids));
                prev_ids = ids;
                long long pc = 0, ec = 0;
                for (const auto& [t, c] : M) pc += is_prefix(ids, t) ? c : 0, ec += t == ids ? c : 0;
                TRIE_CHECK(TR.pass_count(v) == pc && TR.end_count(v) == ec);
                if (v > 0) TRIE_CHECK(TR.child_by_id(TR.parent(v), TR.label(v)) == v);
            }
        }
        // path
        const auto p = TR.path(S);
        TRIE_CHECK(!p.empty() && p[0] == 0);
        for (std::size_t k = 1; k < p.size(); ++k) TRIE_CHECK(TR.parent(p[k]) == p[k - 1] && TR.label(p[k]) == s[k - 1]);
        TRIE_CHECK(p.size() == s.size() + 1 || TR.find(S) == -1);
        ++cases;
    }
}

static void random_test() {
    std::mt19937_64 rng(20261012);
    int cases = 0;
    for (int it = 0; it < 60 * TEST_SCALE; ++it) {
        const int t = it % 4;
        if (t == 0) {
            const Alphabet A = Alphabet::lower();
            const int sg = 1 + (int)(rng() % 3);
            run_ops(rng, A, sg, [&](const Ids& s) { return from_ids(s, A); }, cases);
        } else if (t == 1) {
            const Alphabet A = Alphabet::of("zyx");  // 並べた順が辞書順
            run_ops(rng, A, 3, [&](const Ids& s) { return from_ids(s, A); }, cases);
        } else if (t == 2) {
            const Alphabet A = Alphabet::ints(-5, -2);  // 数の列
            run_ops(rng, A, 3, [&](const Ids& s) {
                std::vector<long long> v;
                for (int c : s) v.push_back(c - 5);
                return v;
            }, cases);
        } else {
            const Alphabet A = Alphabet::ints(2);
            run_ops(rng, A, 2, [&](const Ids& s) { return s; }, cases);
        }
    }
    // 文字列で返すもの
    {
        Trie TR;
        for (const std::string s : {"ab", "a", "abc", "b", "ab"}) TR.insert(s);
        TRIE_CHECK((TR.words() == std::vector<std::string>{"a", "ab", "ab", "abc", "b"}) && TR.kth(3) == "abc" && TR.string_of(TR.find("abc")) == "abc");
        TRIE_CHECK(TR.child(TR.root(), 'a') == TR.find("a") && TR.child(TR.root(), 'z') == -1);
        TR.clear();
        TRIE_CHECK(TR.size() == 0 && TR.node_count() == 1);
        TR.insert("a", 3100000000LL);  // 数が大きくても C(数, 2) を溢れずに
        TRIE_CHECK(TR.sum_lcp_pairs() == 4804999998450000000LL);
        ++cases;
    }
    std::cout << "[trie] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[trie] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    auto ms = [](std::chrono::steady_clock::time_point a) {
        return (long long)std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count();
    };
    std::vector<std::string> vs;
    long long total = 0;
    while (total < 1000000) {
        vs.push_back(rnd_str(rng, 1 + (int)(rng() % 20), 26));
        total += (long long)vs.back().size();
    }
    auto t0 = std::chrono::steady_clock::now();
    Trie TR;
    for (const auto& s : vs) TR.insert(s);
    const long long t_ins = ms(t0);
    t0 = std::chrono::steady_clock::now();
    long long x = 0;
    for (const auto& s : vs) x += TR.count(s) + TR.count_prefix(s.substr(0, 3)) + TR.max_lcp(s, true);
    const long long t_q = ms(t0);
    t0 = std::chrono::steady_clock::now();
    for (int k = 0; k < 100000; ++k) x += (long long)TR.kth(rng() % TR.size()).size() + TR.rank(vs[rng() % vs.size()]);
    const long long t_k = ms(t0);
    t0 = std::chrono::steady_clock::now();
    x += TR.sum_lcp_pairs();
    const long long t_p = ms(t0);
    std::cout << "[trie] benchmark: 文字の合計 10^6（" << vs.size() << " 個）insert " << t_ins << " ms, count + count_prefix + max_lcp " << t_q
              << " ms, kth + rank ×10^5 " << t_k << " ms, sum_lcp_pairs " << t_p << " ms, ノード " << TR.node_count() << "  (check " << x % 7 << ")\n";
}

}  // namespace trie_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace trie_demo {

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
#define SHOW(expr, comment) trie_demo::show(#expr, (expr), comment)
#define RUN(...)                   \
    trie_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("Trie（trie） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  insert / erase / erase_prefix / clear                      追加（cnt 個）・削除・接頭辞で全部消す\n"
              << "  count / contains / count_prefix / count_prefixes_of / size   数・接頭辞の数・接頭辞である数\n"
              << "  kth / kth_ids / rank / words / words_ids                     辞書順 k 番目・順位・全部\n"
              << "  max_lcp / sum_lcp_with / sum_lcp_pairs                       LCP の最大・和\n"
              << "  find / path / child / child_by_id / parent / depth / label / pass_count / end_count / string_of / preorder   ノード\n"
              << "  Trie TR(A) の A で文字の種類（既定は英小文字。Alphabet::upper()・of(\"01\")・ints(K) など。番号の順が辞書順）\n";

    RUN(Trie TR);
    RUN(for (const std::string s : {"apple", "app", "apply", "banana", "app"}) TR.insert(s));

    api("【1】 insert(s, cnt) / erase(s, cnt) / erase_prefix(p)", "int（ノード）/ bool / long long",
        "cnt 個入れる（0 個でもノードは作る）/ cnt 個消す（足りなければ何もせず false）/ p で始まる文字列を全部消す（消した数）", "O(|s|)");
    SHOW(TR.size(), "app は 2 個");
    SHOW(TR.erase("ap"), "入っていない → false");

    api("【2】 count(s) / contains(s) / count_prefix(p) / count_prefixes_of(s)", "long long / bool",
        "s の数 / 入っているか / p で始まる文字列の数（p 自身も）/ s の接頭辞である文字列の数（s 自身も）", "O(|s|)");
    SHOW(TR.count("app"), "");
    SHOW(TR.contains("appl"), "接頭辞だけでは入っていない");
    SHOW(TR.count_prefix("app"), "app, app, apple, apply");
    SHOW(TR.count_prefixes_of("apples"), "app, app, apple");

    api("【3】 kth(k) / rank(s) / words()", "string / long long / vector<string>", "辞書順 k 番目（0 から。同じものも数える）/ s より小さい数 / 全部を辞書順に", "O(|答え| σ) / O(|s| σ)");
    SHOW(TR.kth(2), "app, app, apple, ...");
    SHOW(TR.rank("apply"), "");
    SHOW(TR.words(), "");

    api("【4】 max_lcp(s, exclude_self) / sum_lcp_with(s) / sum_lcp_pairs()", "int / long long / long long",
        "入っている文字列との LCP の最大（exclude_self = true で s を 1 つ除く）/ 全部との LCP の和 / 全部の組（i < j）の LCP の和", "O(|s|) / O(ノード数)");
    SHOW(TR.max_lcp("applause"), "appl");
    SHOW(TR.max_lcp("banana", true), "自分を除くと共通なし");
    SHOW(TR.sum_lcp_with("apx"), "2 + 2 + 2 + 2 + 0");
    SHOW(TR.sum_lcp_pairs(), "");

    api("【5】 ノード：find(s) / path(s) / child(v, c) / parent(v) / depth(v) / label(v) / pass_count(v) / end_count(v) / string_of(v) / preorder()",
        "int / vector<int> / ...", "ノード 0 が根。find はなければ -1（数が 0 でもノードがあれば返す）。path は根から s をたどったノードの列。preorder は辞書順の行きがけ順（逆にすれば子が先。DP 用）", "O(|s|) / O(1)");
    {
        RUN(const int v = TR.find("app"));
        SHOW(TR.pass_count(v), "");
        SHOW(TR.end_count(v), "");
        SHOW(TR.depth(v), "");
        SHOW(TR.string_of(TR.child(v, 'l')), "");
        SHOW(TR.path("apex"), "a, ap まで（e がない）");
    }
    SHOW(TR.erase_prefix("app"), "app で始まる 4 個を消す");
    SHOW(TR.words(), "");

    head("問題での使い方");
    {
        Trie K;
        const std::vector<std::string> S = {"abc", "abb", "aac"};
        for (const auto& s : S) K.insert(s);
        std::vector<int> ans;
        for (const auto& s : S) ans.push_back(K.max_lcp(s, true));
        note("ABC287 E（Karuta）: 全部 insert して max_lcp(S_i, true)。入力例 1 → " + join(ans));
    }
    {
        Trie K;
        for (const std::string s : {"ab", "abc", "arc"}) K.insert(s);
        note("ABC353 E（Yet Another Sigma Problem）: 全部 insert して sum_lcp_pairs()。入力例 1 → " + to_s(K.sum_lcp_pairs()));
    }
    {
        note("ABC268 G（Random Student ID）: 自分の接頭辞の数 P = count_prefixes_of(S_i) - 1、自分で始まる数 Q = count_prefix(S_i) - 1 として");
        note("  期待値 = 1 + P + (N - 1 - P - Q) / 2（mod は math_v2 の ModInt998 で）");
        Trie K;
        const std::vector<std::string> S = {"a", "aa", "ab"};
        for (const auto& s : S) K.insert(s);
        std::vector<std::string> ans;
        for (const auto& s : S) {
            const long long P = K.count_prefixes_of(s) - 1, Q = K.count_prefix(s) - 1, N = (long long)S.size();
            const long long twice = 2 + 2 * P + (N - 1 - P - Q);  // 2 倍した値
            ans.push_back(twice % 2 ? to_s(twice) + "/2" : to_s(twice / 2));
        }
        note("  入力例 1 → " + join(ans) + "（5/2 は mod 998244353 で 499122179）");
    }
    {
        note("ABC377 G（Edit to Match）: ノードごとに best[v] = min(|S_j| - depth(v)) を持ち、S_k の path の上で |S_k| - depth(v) + best[v] の最小");
        Trie K;
        std::vector<long long> best;
        std::vector<long long> ans;
        for (const std::string s : {"snuke", "snuki", "snuuk"}) {
            const long long n = (long long)s.size();
            long long a = n;
            for (int v : K.path(s))
                if (v < (int)best.size()) a = std::min(a, n - K.depth(v) + best[v]);
            ans.push_back(a);
            K.insert(s);
            best.resize(K.node_count(), (long long)4e18);
            for (int v : K.path(s)) best[v] = std::min(best[v], n - K.depth(v));
        }
        note("  入力例 1 → " + join(ans));
    }
    {
        note("ABC403 E（Forbidden Prefix）: X は insert(X, 0) でノードを作って印を付け、erase_prefix(X)。Y は path に印がなければ insert。答えは size()");
        Trie K;
        std::vector<char> banned;
        std::vector<long long> ans;
        const std::vector<std::pair<int, std::string>> qs = {{1, "at"}, {2, "watcoder"}, {2, "atcoder"}, {1, "wa"}};
        for (const auto& [t, s] : qs) {
            if (t == 1) {
                const int v = K.insert(s, 0);
                banned.resize(K.node_count(), 0);
                banned[v] = 1;
                K.erase_prefix(s);
            } else {
                bool ok = true;
                for (int v : K.path(s)) ok = ok && !(v < (int)banned.size() && banned[v]);
                if (ok) K.insert(s);
            }
            ans.push_back(K.size());
        }
        note("  入力例 1 → " + join(ans));
    }
}

}  // namespace trie_demo

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
    trie_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    trie_test::benchmark();    // 速度計測
    trie_demo::demo();         // 使い方デモ
    return 0;
}
