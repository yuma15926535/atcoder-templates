// =============================================================================
//  ローリングハッシュ（部分文字列の一致・LCP・辞書順の比較・連結・2 次元）   --- 文字列第二世代
// -----------------------------------------------------------------------------
//  分類       : 01 ハッシュ
//  依存       : 00_str_base（同梱済み）
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）。unsigned __int128 を使う（GCC / Clang）
//  文字コード : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  方式：mod 2^61 - 1 の 1 本。基数は実行ごとの乱数で、すべての文字列で共通（別の文字列の部分文字列どうし、
//        つなげたものどうしも比べられる）。累乗の表は基数を変えずに伸ばす。1 回の比較で衝突する確率は約 N / 2^61
//  列は string でも vector<int>・vector<long long> でもよい（値はそのまま区別する。10^18 級・負の値も）
//  範囲は半開区間 [l, r)。「L 文字目から R 文字目」は (L - 1, R)
//
//  機能一覧                                                                       計算量
//    StrHash{value, len}                   部分文字列のハッシュ（長さも持つ）
//      a + b / a += b                      つなげたもの（途中に挿入・削除した文字列も作れる）          O(1)
//      a.repeat(k)                         k 回つなげたもの                                          O(log k)
//      == / != / <（map・sort 用）          unordered_map<StrHash, V> もそのまま使える
//    hash_of(s) / hash_of_char(c) / empty_hash()   列全体 / 1 文字 / 空のハッシュ                       O(N) / O(1)
//    RollingHash H(s)                      前処理 O(N)
//      H.get(l, r) / get_len(l, len)       s[l, r) / s[l, l + len) のハッシュ                       O(1)
//      H.prefix(r) / suffix(l) / all()     s[0, r) / s[l, N) / 全体
//      H.same(l1, r1, l2, r2)              s[l1, r1) == s[l2, r2) か                                 O(1)
//      H.lcp(i, j) / H.lcp(i, T, j)        s[i..) と s[j..)（t[j..)）の最長共通接頭辞の長さ          O(log N)
//      H.common_suffix(i, j) / (i, T, j)   s[0, i) と s[0, j)（t[0, j)）の最長共通接尾辞の長さ        O(log N)
//      H.compare(l1, r1, l2, r2) / (l1, r1, T, l2, r2)   辞書順の比較（-1 / 0 / 1）                  O(log N)
//      H.find_all(pat) / count(pat)        pat の出現位置（重なりも）/ 数（pat は列でも StrHash でも）   O(N)
//      H.code(i) / size()                  i 文字目の値（char は 0〜255）/ 長さ
//    RollingHash2D G(grid)                 グリッドの長方形のハッシュ（vector<string> でも vector<vector<int>> でも）  前処理 O(HW)
//      G.get(i1, j1, i2, j2)               行 [i1, i2) × 列 [j1, j2) のハッシュ（GridHash{value, h, w}）  O(1)
//      G.find_all(pat)                     pat と一致する左上の位置 (i, j) の一覧                     O(HW)
//      G.all() / rows() / cols()
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    rh_test::random_test()  … 愚直な比較（std::string の比較・全部の位置）と全機能を突き合わせるランダムテスト
//    rh_test::benchmark()    … 速度計測
//    rh_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
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
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト
// =============================================================================
namespace rh_test {

#define RH_CHECK(...)                                                                          \
    do {                                                                                        \
        if (!(__VA_ARGS__)) {                                                                   \
            std::cerr << "[rolling_hash] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
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

static int bf_lcp(const std::string& a, int i, const std::string& b, int j) {
    int k = 0;
    while (i + k < (int)a.size() && j + k < (int)b.size() && a[i + k] == b[j + k]) ++k;
    return k;
}
static int bf_cmp(const std::string& a, const std::string& b) { return a < b ? -1 : (a == b ? 0 : 1); }

static void test_1d(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 400 * TEST_SCALE; ++it) {
        const int sg = 1 + (int)(rng() % 3);
        const std::string s = rnd_str_shaped(rng, (int)(rng() % 14), sg), t = rnd_str_shaped(rng, (int)(rng() % 14), sg);
        const RollingHash H(s), T(t);
        const int n = (int)s.size(), m = (int)t.size();
        RH_CHECK(H.size() == n && H.all() == hash_of(s) && T.all() == hash_of(t));
        for (int l1 = 0; l1 <= n; ++l1)
            for (int r1 = l1; r1 <= n; ++r1) {
                const std::string a = s.substr(l1, r1 - l1);
                RH_CHECK(H.get(l1, r1) == hash_of(a) && H.get_len(l1, r1 - l1) == H.get(l1, r1));
                for (int mid = l1; mid <= r1; ++mid) RH_CHECK(H.get(l1, mid) + H.get(mid, r1) == H.get(l1, r1));
                for (int l2 = 0; l2 <= m; ++l2)
                    for (int r2 = l2; r2 <= m; ++r2) {
                        const std::string b = t.substr(l2, r2 - l2);
                        RH_CHECK((H.get(l1, r1) == T.get(l2, r2)) == (a == b));
                        RH_CHECK(H.compare(l1, r1, T, l2, r2) == bf_cmp(a, b));
                        ++cases;
                    }
                for (int l2 = 0; l2 <= n; ++l2) {
                    const int r2 = l2 + (int)(rng() % (n - l2 + 1));
                    const std::string b = s.substr(l2, r2 - l2);
                    RH_CHECK(H.same(l1, r1, l2, r2) == (a == b) && H.compare(l1, r1, l2, r2) == bf_cmp(a, b));
                }
            }
        for (int i = 0; i <= n; ++i) {
            for (int j = 0; j <= n; ++j) {
                RH_CHECK(H.lcp(i, j) == bf_lcp(s, i, s, j));
                int k = 0;
                while (k < i && k < j && s[i - 1 - k] == s[j - 1 - k]) ++k;
                RH_CHECK(H.common_suffix(i, j) == k);
            }
            for (int j = 0; j <= m; ++j) {
                RH_CHECK(H.lcp(i, T, j) == bf_lcp(s, i, t, j));
                int k = 0;
                while (k < i && k < j && s[i - 1 - k] == t[j - 1 - k]) ++k;
                RH_CHECK(H.common_suffix(i, T, j) == k);
            }
        }
        // 出現位置
        std::vector<int> ex;
        for (int i = 0; i + m <= n; ++i)
            if (s.compare(i, m, t) == 0) ex.push_back(i);
        RH_CHECK(H.find_all(t) == ex && H.find_all(hash_of(t)) == ex && H.count(t) == (long long)ex.size());
        // 連結・繰り返し・1 文字
        RH_CHECK(hash_of(s) + hash_of(t) == hash_of(s + t) && hash_of(s + t) == hash_of(s) + T.all());
        StrHash acc = empty_hash();
        for (char c : s) acc += hash_of_char(c);
        RH_CHECK(acc == hash_of(s));
        const long long k = rng() % 6;
        RH_CHECK(hash_of(s).repeat(k) == hash_of(repeat_string(s, k)));
        RH_CHECK(hash_of(s).repeat(0) == empty_hash());
        // 途中に挿入した文字列のハッシュ
        const int p = (int)(rng() % (n + 1));
        RH_CHECK(H.prefix(p) + T.all() + H.suffix(p) == hash_of(s.substr(0, p) + t + s.substr(p)));
    }
    // 巨大な繰り返し（長さだけでなく値も：(a^k)^2 = a^(2k)）
    const StrHash a = hash_of("abc");
    const StrHash big = a.repeat(1000000000000LL);
    RH_CHECK(big.len == 3000000000000LL && big + big == a.repeat(2000000000000LL) && big.repeat(3) == a.repeat(3000000000000LL));
}

static void test_values(std::mt19937_64& rng, int& cases) {
    // 数の列：大きな値・負の値・2^61 - 1 の倍数ずれも区別する
    const long long M = (long long)((1ULL << 61) - 1);
    const std::vector<long long> sp = {0, 1, -1, M, M + 1, 2 * M, -M, (1LL << 60) - 1, 1LL << 60, std::numeric_limits<long long>::max(), std::numeric_limits<long long>::min()};
    for (int it = 0; it < 2000 * TEST_SCALE; ++it) {
        const int n = 1 + (int)(rng() % 6);
        std::vector<long long> a(n), b(n);
        for (auto& x : a) x = rng() % 2 ? sp[rng() % sp.size()] : rnd_ll(rng);
        b = a;
        if (rng() % 2) b[rng() % n] = rng() % 2 ? sp[rng() % sp.size()] : rnd_ll(rng);
        RH_CHECK((hash_of(a) == hash_of(b)) == (a == b));
        const RollingHash A(a), B(b);
        int k = 0;
        while (k < n && a[k] == b[k]) ++k;
        RH_CHECK(A.lcp(0, B, 0) == k);
        RH_CHECK(A.compare(0, n, B, 0, n) == (a < b ? -1 : (a == b ? 0 : 1)));
        ++cases;
    }
    // string と vector<int> は同じ値なら同じハッシュ
    const std::string s = "hello";
    RH_CHECK(hash_of(s) == hash_of(std::vector<int>(s.begin(), s.end())) && hash_of("hello") == hash_of(s));
    // 128 以上の char
    std::string hi = "a";
    hi += (char)200;
    RH_CHECK(RollingHash(hi).compare(0, 1, 1, 2) == -1);
}

static void test_distinct(std::mt19937_64& rng, int& cases) {
    // 異なる部分文字列の数をハッシュと set<string> で（衝突しないことの確認）
    for (int it = 0; it < 6 * TEST_SCALE; ++it) {
        const std::string s = rnd_str_shaped(rng, 250, 2 + (int)(rng() % 2));
        const RollingHash H(s);
        std::unordered_set<StrHash> hs;
        std::set<std::string> ss;
        for (int l = 0; l < (int)s.size(); ++l)
            for (int r = l + 1; r <= (int)s.size(); ++r) hs.insert(H.get(l, r)), ss.insert(s.substr(l, r - l));
        RH_CHECK(hs.size() == ss.size());
        cases += (int)ss.size();
    }
    std::unordered_map<StrHash, int> mp;
    mp[hash_of("ab")] = 3;
    RH_CHECK(mp[hash_of("a") + hash_of("b")] == 3);
}

static void test_2d(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 300 * TEST_SCALE; ++it) {
        const int H = (int)(rng() % 6), W = H ? 1 + (int)(rng() % 6) : 0, sg = 1 + (int)(rng() % 2);
        std::vector<std::string> g(H);
        for (auto& row : g) row = rnd_str(rng, W, sg);
        const RollingHash2D G(g);
        RH_CHECK(G.rows() == H && G.cols() == W);
        auto sub = [&](int i1, int j1, int i2, int j2) {
            std::vector<std::string> r;
            for (int i = i1; i < i2; ++i) r.push_back(g[i].substr(j1, j2 - j1));
            return r;
        };
        for (int t = 0; t < 30; ++t) {
            const int i1 = (int)(rng() % (H + 1)), i2 = i1 + (int)(rng() % (H - i1 + 1));
            const int j1 = (int)(rng() % (W + 1)), j2 = j1 + (int)(rng() % (W - j1 + 1));
            const int k1 = (int)(rng() % (H - (i2 - i1) + 1)), l1 = (int)(rng() % (W - (j2 - j1) + 1));
            const auto a = sub(i1, j1, i2, j2), b = sub(k1, l1, k1 + (i2 - i1), l1 + (j2 - j1));
            RH_CHECK((G.get(i1, j1, i2, j2) == G.get(k1, l1, k1 + i2 - i1, l1 + j2 - j1)) == (a == b || i1 == i2 || j1 == j2));
            ++cases;
        }
        // 模様探し
        const int h = (int)(rng() % 3) + 1, w = (int)(rng() % 3) + 1;
        std::vector<std::string> pat(h);
        for (auto& row : pat) row = rnd_str(rng, w, sg);
        std::vector<std::pair<int, int>> ex;
        for (int i = 0; i + h <= H; ++i)
            for (int j = 0; j + w <= W; ++j)
                if (sub(i, j, i + h, j + w) == pat) ex.push_back({i, j});
        RH_CHECK(G.find_all(pat) == ex);
        std::vector<std::vector<int>> gi(H, std::vector<int>(W));
        for (int i = 0; i < H; ++i)
            for (int j = 0; j < W; ++j) gi[i][j] = g[i][j];
        RH_CHECK(RollingHash2D(gi).all() == G.all());
    }
}

static void test_changes(std::mt19937_64& rng, int& cases) {  // center_palindrome_lengths_with_changes
    for (int it = 0; it < 1500 * TEST_SCALE; ++it) {
        const int n = (int)(rng() % 25), k = (int)(rng() % 4);
        const std::string s = rnd_str_shaped(rng, n, 1 + (int)(rng() % 3));
        const auto L = center_palindrome_lengths_with_changes(s, k);
        RH_CHECK((int)L.size() == std::max(0, 2 * n - 1));
        for (int c = 0; c < 2 * n - 1; ++c) {
            int best = 0;
            for (int len = c % 2 == 0 ? 1 : 2; len <= n; len += 2) {  // 中心 c の長さ len の区間
                const int l = (c + 1 - len) / 2;
                if (l < 0 || l + len > n) break;
                int mis = 0;
                for (int t = 0; t < len / 2; ++t) mis += s[l + t] != s[l + len - 1 - t];
                if (mis <= k) best = len;
                else break;
            }
            RH_CHECK(L[c] == best);
        }
        if (k == 0) {  // 04_palindrome の center_palindrome_lengths と同じ（愚直に）
            for (int c = 0; c < 2 * n - 1; ++c) {
                int len = c % 2 == 0 ? 1 : 0, l = c % 2 == 0 ? c / 2 - 1 : (c - 1) / 2, r = c % 2 == 0 ? c / 2 + 1 : (c + 1) / 2;
                while (l >= 0 && r < n && s[l] == s[r]) --l, ++r, len += 2;
                RH_CHECK(L[c] == len);
            }
        }
        ++cases;
    }
}

static void random_test() {
    std::mt19937_64 rng(20261012);
    int cases = 0;
    test_1d(rng, cases);
    test_values(rng, cases);
    test_distinct(rng, cases);
    test_2d(rng, cases);
    test_changes(rng, cases);
    std::cout << "[rolling_hash] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[rolling_hash] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    auto ms = [](std::chrono::steady_clock::time_point a) {
        return (long long)std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count();
    };
    const int n = 10000000;
    const std::string s = rnd_str(rng, n, 26);
    auto t0 = std::chrono::steady_clock::now();
    const RollingHash H(s);
    const long long t_build = ms(t0);
    std::vector<int> ls(n), rs(n);
    for (int i = 0; i < n; ++i) {
        ls[i] = (int)(rng() % (n + 1)), rs[i] = (int)(rng() % (n + 1));
        if (ls[i] > rs[i]) std::swap(ls[i], rs[i]);
    }
    t0 = std::chrono::steady_clock::now();
    unsigned long long x = 0;
    for (int i = 0; i < n; ++i) x += H.get(ls[i], rs[i]).value;
    const long long t_get = ms(t0);
    const std::string s2 = rnd_str_shaped(rng, 1000000, 2);
    const RollingHash H2(s2);
    t0 = std::chrono::steady_clock::now();
    long long y = 0;
    for (int i = 0; i < 1000000; ++i) y += H2.lcp((int)(rng() % 1000001), (int)(rng() % 1000001));
    const long long t_lcp = ms(t0);
    t0 = std::chrono::steady_clock::now();
    const auto occ = H.find_all(s.substr(5000000, 3));
    const long long t_find = ms(t0);
    std::vector<std::string> g(2000, std::string(2000, 'a'));
    for (auto& row : g)
        for (char& c : row) c = (char)('a' + rng() % 2);
    t0 = std::chrono::steady_clock::now();
    const RollingHash2D G(g);
    const auto f2 = G.find_all(std::vector<std::string>{"ab", "ba"});
    const long long t_2d = ms(t0);
    std::cout << "[rolling_hash] benchmark: 構築（N = 10^7）" << t_build << " ms, get × 10^7 " << t_get << " ms, lcp × 10^6（N = 10^6）" << t_lcp
              << " ms, find_all（N = 10^7）" << t_find << " ms, 2D 構築 + find_all（2000×2000）" << t_2d << " ms  (check "
              << (x + (unsigned long long)y + occ.size() + f2.size()) % 7 << ")\n";
}

}  // namespace rh_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace rh_demo {

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
#define SHOW(expr, comment) rh_demo::show(#expr, (expr), comment)
#define RUN(...)                   \
    rh_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("ローリングハッシュ（rolling_hash） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  RollingHash H(s)：get / get_len / prefix / suffix / all / same / lcp / common_suffix / compare / find_all / count\n"
              << "  StrHash{value, len}：+（連結）/ repeat(k) / == / unordered_map のキーに\n"
              << "  hash_of(s) / hash_of_char(c) / empty_hash()\n"
              << "  RollingHash2D G(grid)：get(i1, j1, i2, j2) / find_all(pat) / all\n"
              << "  center_palindrome_lengths_with_changes(s, k)   k 文字まで書き換えてよい回文（中心ごとの最長）\n"
              << "  基数は実行ごとの乱数で全部の文字列で共通なので、別の文字列の部分文字列どうしも比べられる\n";

    api("【1】 RollingHash H(s) / H.get(l, r) / H.same(l1, r1, l2, r2)", "StrHash / bool",
        "s[l, r) のハッシュ。同じ長さで値が同じなら（ほぼ確実に）同じ文字列。何度でも O(1)", "構築 O(N)・取得 O(1)");
    RUN(const std::string S = "abcabcab");
    RUN(const RollingHash H(S));
    SHOW(H.same(0, 3, 3, 6), "\"abc\" と \"abc\"");
    SHOW(H.same(0, 3, 5, 8), "\"abc\" と \"cab\"");
    SHOW(H.get(1, 4) == hash_of("bca"), "別に作ったハッシュとも比べられる");
    note("問題の「L 文字目から R 文字目」は H.get(L - 1, R)。");

    api("【2】 H.lcp(i, j) / H.common_suffix(i, j) / H.compare(l1, r1, l2, r2)", "int",
        "s[i..) と s[j..) の最長共通接頭辞・s[0, i) と s[0, j) の最長共通接尾辞・部分文字列の辞書順（-1 / 0 / 1）。二分探索", "O(log N)");
    SHOW(H.lcp(0, 3), "\"abcabcab\" と \"abcab\" → 5");
    SHOW(H.common_suffix(8, 5), "\"abcabcab\" と \"abcab\" の共通の末尾");
    SHOW(H.compare(1, 3, 0, 3), "\"bc\" と \"abc\" → 1（大きい）");
    RUN(const RollingHash T("abd"));
    SHOW(H.lcp(0, T, 0), "別の文字列 \"abd\" とも");
    SHOW(H.compare(0, 3, T, 0, 3), "\"abc\" < \"abd\"");
    note("O(1) の LCP・比較が要るなら 03_suffix_array の SuffixArray。");

    api("【3】 StrHash の +（連結）/ repeat(k) / hash_of(s) / hash_of_char(c)", "StrHash",
        "ハッシュどうしをつなげる（文字列を作らない）。途中に挿入・削除・入れ替えた文字列と比べるときに", "O(1) / O(log k)");
    SHOW(H.get(0, 2) + H.get(5, 8) == hash_of("abcab"), "\"ab\" + \"cab\"");
    SHOW(hash_of("ab").repeat(3) == hash_of("ababab"), "3 回つなげる");
    SHOW(hash_of("x").repeat(1000000000000LL).len, "長さ 10^12 の文字列のハッシュも O(log)");
    {
        // ABC284 F: T（長さ 2N）= S[0, i) + rev(S) + S[i, N) となる S と i
        const std::string Tt = "abcbac";
        const int N = 3;
        const RollingHash F(Tt), R(reversed(Tt));
        int ans = -1;
        for (int i = 0; i <= N && ans < 0; ++i)  // T[0, i) + T[i + N, 2N) と、T[i, i + N) を反転したもの
            if (F.prefix(i) + F.suffix(i + N) == R.get(N - i, 2 * N - i)) ans = i;
        note("ABC284 F（ABCBAC）: F.prefix(i) + F.suffix(i + N) と、反転した T のハッシュ R.get(N - i, 2N - i) を比べる");
        std::cout << "    T = abcbac → i = " << ans << "（S = " << Tt.substr(0, ans) + Tt.substr(ans + N) << "）\n";
    }

    api("【4】 unordered_map<StrHash, V> / H.find_all(pat) / H.count(pat)", "—",
        "ハッシュを map のキーに（同じ部分文字列をまとめる）。find_all は pat の出現位置（重なりも）", "O(1) / O(N)");
    {
        // ABC141 E: 重ならずに 2 回以上現れる部分文字列の最長（長さで二分探索）
        const std::string S2 = "ababa";
        const RollingHash G(S2);
        const int n = (int)S2.size();
        auto ok = [&](int len) {
            std::unordered_map<StrHash, int> first;  // ハッシュ → 最初の開始位置
            for (int i = 0; i + len <= n; ++i) {
                auto [it, ins] = first.emplace(G.get(i, i + len), i);
                if (!ins && i - it->second >= len) return true;
            }
            return false;
        };
        int lo = 0, hi = n / 2;
        while (lo < hi) {
            const int mid = (lo + hi + 1) / 2;
            if (ok(mid)) lo = mid;
            else hi = mid - 1;
        }
        note("ABC141 E（Who Says a Pun?）: 長さ len で二分探索し、unordered_map<StrHash, int> に最初の位置を覚える");
        std::cout << "    S = ababa → " << lo << "\n";
    }
    SHOW(H.find_all("ab"), "\"ab\" の出現位置");
    SHOW(H.count("abcab"), "重なりも数える");

    api("【5】 数の列にも使える", "—", "vector<int>・vector<long long> の部分列の一致（10^18 級・負の値も区別する）", "—");
    RUN(const std::vector<long long> A = {5, -1, 1000000000000000000LL, 5, -1});
    RUN(const RollingHash HA(A));
    SHOW(HA.same(0, 2, 3, 5), "{5, -1} と {5, -1}");

    api("【6】 RollingHash2D G(grid) / G.get(i1, j1, i2, j2) / G.find_all(pat)", "GridHash / vector<pair<int,int>>",
        "グリッドの長方形 [i1, i2) × [j1, j2) のハッシュ。大きなグリッドの中で模様を探す（小さいなら graph_grid_v2 の find_subgrid）", "構築 O(HW)・取得 O(1)");
    RUN(const std::vector<std::string> g = {"#..#", ".##.", "#..#"});
    RUN(const RollingHash2D G(g));
    SHOW(G.get(0, 0, 1, 2) == G.get(2, 2, 3, 4), "\"#.\" と \".#\" は違う");
    SHOW(G.get(0, 0, 1, 4) == G.get(2, 0, 3, 4), "1 行目と 3 行目は同じ");
    SHOW(G.find_all(std::vector<std::string>{"#.", ".#"}).size(), "2×2 の模様の数");
    SHOW(G.find_all(std::vector<std::string>{"#.", ".#"})[0].second, "最初の位置の列");

    api("【7】 center_palindrome_lengths_with_changes(s, k)", "vector<int>（長さ 2N-1）",
        "中心 c ごとに、食い違う組が k 個以下の最長の長さ（偶数 c は s[c/2] が中心の奇数長、奇数 c は間が中心の偶数長）。中心から LCP で k + 1 回広げる", "O(N k log N)");
    SHOW(center_palindrome_lengths_with_changes("abcda", 1), "1 文字まで書き換えてよい");
    {
        // ABC468 D 型：1 文字まで書き換えて回文にできる部分文字列の数 = 中心ごとに (L + 1) / 2 を足す（同じ偶奇の長さは全部できる）
        const auto L = center_palindrome_lengths_with_changes(std::string("abca"), 1);
        long long cnt = 0;
        for (int c = 0; c < (int)L.size(); ++c) cnt += c % 2 == 0 ? (L[c] + 1) / 2 : L[c] / 2;
        note("1 文字まで書き換えて回文になる部分文字列の数：中心ごとに長さ L 以下で同じ偶奇のもの（奇数中心は (L+1)/2 個、偶数中心は L/2 個）");
        std::cout << "    S = abca → " << cnt << "（a, b, c, a, ab, bc, ca, abc, bca, abca の 10 個）\n";
    }
}

}  // namespace rh_demo

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
    rh_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    rh_test::benchmark();    // 速度計測
    rh_demo::demo();         // 使い方デモ
    return 0;
}
