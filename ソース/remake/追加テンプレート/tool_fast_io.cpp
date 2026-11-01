// =============================================================================
//  高速入出力（fread でまとめて読む・fwrite でまとめて書く。入力が 10^6 個を超えるとき）   --- AtCoder 用テンプレート
// -----------------------------------------------------------------------------
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）
//  文字コード     : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  cin（sync_with_stdio(false) にしても）より数倍速く、整数・文字列・実数を読み書きする。
//  入力の数が 10^6〜10^7 個のときや、出力の行がとても多いときに使う
//
//  機能一覧                                                                                  計算量
//    fast_in.read<T>() / fast_in.read(x, y, …)   1 つ読む / まとめて読む。T は整数（int・long long・unsigned・__int128）・
//                                                char（空白以外の 1 文字）・string（空白まで）・double / long double     O(文字数)
//    fast_in >> x >> y                           cin と同じ書き方
//    fast_in.read_vec<T>(n)                      n 個読んで vector<T>
//    fast_in.eof()                               もう読むものが無いか（空白を飛ばしてから判定）
//    fast_out.write(x)                           1 つ書く（区切りなし）。整数・char・文字列・double（桁数は set_precision）
//    fast_out.print(a, b, …)                     空白区切りで書いて改行（vector は中身を空白区切り）
//    fast_out << x                               cout と同じ書き方（'\n' も書ける）
//    fast_out.set_precision(k) / fast_out.flush()   小数点以下の桁数（既定 15）/ たまっている分を書き出す
//    FastInput in(fp) / FastOutput out(fp)       ファイル（FILE*）を指定して作る（fast_in / fast_out は stdin / stdout）
//
//  注意
//    fast_in と cin、fast_out と cout を混ぜない（読む・書く順番がずれる）。どうしても混ぜるなら fast_out.flush() してから cout
//    fast_out はプログラムの終わり（main から return・exit）で自動で書き出す。abort・quick_exit では消えるので flush する
//    インタラクティブ問題では使わない（tool_interactive.cpp を使う）
//  速さの目安: 10^7 個の整数を読む 約 0.17 秒（istream の >> だと 約 1.3 秒）、書く 約 0.3 秒（一時ファイルに）
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    fastio_test::random_test()  … 一時ファイルに書いた入力を stringstream の結果と突き合わせ
//    fastio_test::benchmark()    … 速度計測（cin・cout とも比べる）
//    fastio_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
// =============================================================================
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <random>
#include <sstream>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

// ======================= ここからテンプレート本体（コピー範囲） =======================
namespace fastio_internal {
template <class T>
struct is_vec : std::false_type {};
template <class T, class A>
struct is_vec<std::vector<T, A>> : std::true_type {};
}  // namespace fastio_internal

struct FastInput {
    explicit FastInput(FILE* fp = stdin) : fp_(fp) {}
    template <class T>
    T read() {
        T x{};
        get(x);
        return x;
    }
    template <class... Args>
    void read(Args&... args) {
        (get(args), ...);
    }
    template <class T>
    std::vector<T> read_vec(int n) {
        std::vector<T> v(n);
        for (auto& x : v) get(x);
        return v;
    }
    template <class T>
    FastInput& operator>>(T& x) {
        get(x);
        return *this;
    }
    bool eof() {
        skip();
        return peek() < 0;
    }

   private:
    FILE* fp_;
    static constexpr size_t SZ = 1 << 16;
    char buf_[SZ];
    size_t pos_ = 0, len_ = 0;
    int peek() {
        if (pos_ == len_) {
            len_ = std::fread(buf_, 1, SZ, fp_), pos_ = 0;
            if (len_ == 0) return -1;
        }
        return (unsigned char)buf_[pos_];
    }
    int next() {
        const int c = peek();
        if (c >= 0) ++pos_;
        return c;
    }
    void skip() {
        int c;
        while ((c = peek()) >= 0 && c <= ' ') ++pos_;
    }
    std::string token() {
        skip();
        std::string s;
        int c;
        while ((c = peek()) > ' ') s += (char)c, ++pos_;
        return s;
    }
    template <class T>
    void get(T& x) {
        if constexpr (std::is_same_v<T, char>) {
            skip();
            const int c = next();
            x = c < 0 ? '\0' : (char)c;
        } else if constexpr (std::is_same_v<T, std::string>) {
            x = token();
        } else if constexpr (std::is_floating_point_v<T>) {
            const std::string s = token();
            x = (T)std::strtold(s.c_str(), nullptr);
        } else {  // 整数（__int128 も）
            skip();
            bool neg = false;
            if (peek() == '-') neg = true, ++pos_;
            else if (peek() == '+') ++pos_;
            using U = std::conditional_t<(sizeof(T) > 8), unsigned __int128, unsigned long long>;
            U v = 0;
            int c;
            while ((c = peek()) >= '0' && c <= '9') v = v * 10 + (U)(c - '0'), ++pos_;
            x = neg ? (T)(U(0) - v) : (T)v;
        }
    }
};

struct FastOutput {
    explicit FastOutput(FILE* fp = stdout) : fp_(fp) {}
    ~FastOutput() { flush(); }
    FastOutput(const FastOutput&) = delete;
    FastOutput& operator=(const FastOutput&) = delete;
    void set_precision(int k) { prec_ = k; }
    void flush() {
        if (len_) std::fwrite(buf_, 1, len_, fp_), len_ = 0;
        std::fflush(fp_);
    }
    template <class T>
    void write(const T& x) {
        if constexpr (std::is_same_v<T, char>) {
            put(x);
        } else if constexpr (std::is_same_v<T, std::string>) {
            for (char c : x) put(c);
        } else if constexpr (std::is_same_v<std::decay_t<T>, const char*> || std::is_same_v<std::decay_t<T>, char*>) {
            for (const char* p = x; *p; ++p) put(*p);
        } else if constexpr (std::is_array_v<T>) {  // 文字列リテラル
            for (const char* p = x; *p; ++p) put(*p);
        } else if constexpr (std::is_same_v<T, bool>) {
            put(x ? '1' : '0');
        } else if constexpr (std::is_floating_point_v<T>) {
            char tmp[96];
            const int n = std::snprintf(tmp, sizeof tmp, "%.*Lf", prec_, (long double)x);
            for (int i = 0; i < n && i < (int)sizeof tmp - 1; ++i) put(tmp[i]);
        } else {  // 整数
            using U = std::conditional_t<(sizeof(T) > 8), unsigned __int128, unsigned long long>;
            U v = (U)x;
            if constexpr (std::is_signed_v<T> || std::is_same_v<T, __int128>) {
                if (x < 0) put('-'), v = U(0) - v;
            }
            char tmp[48];
            int n = 0;
            do tmp[n++] = (char)('0' + (int)(v % 10)), v /= 10;
            while (v);
            while (n) put(tmp[--n]);
        }
    }
    template <class... Args>
    void print(const Args&... args) {
        bool first = true;
        (one(args, first), ...);
        (void)first;
        put('\n');
    }
    template <class T>
    FastOutput& operator<<(const T& x) {
        write(x);
        return *this;
    }

   private:
    FILE* fp_;
    static constexpr size_t SZ = 1 << 16;
    char buf_[SZ];
    size_t len_ = 0;
    int prec_ = 15;
    void put(char c) {
        if (len_ == SZ) std::fwrite(buf_, 1, len_, fp_), len_ = 0;
        buf_[len_++] = c;
    }
    template <class T>
    void one(const T& x, bool& first) {
        if constexpr (fastio_internal::is_vec<T>::value) {
            for (const auto& e : x) one(e, first);
        } else {
            if (!first) put(' ');
            write(x);
            first = false;
        }
    }
};
inline FastInput fast_in;
inline FastOutput fast_out;
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト・速度計測
// =============================================================================
namespace fastio_test {

#if defined(_DEBUG) || defined(_GLIBCXX_DEBUG) || (defined(__GNUC__) && !defined(__OPTIMIZE__))
[[maybe_unused]] static const int TEST_SCALE = 1;
#else
[[maybe_unused]] static const int TEST_SCALE = 8;
#endif

#define FI_CHECK(...)                                                                          \
    do {                                                                                       \
        if (!(__VA_ARGS__)) {                                                                  \
            std::cerr << "[tool_fast_io] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
            std::exit(1);                                                                      \
        }                                                                                      \
    } while (0)

// [lo, hi] の一様乱数
[[maybe_unused]] static long long rnd(std::mt19937_64& rng, long long lo, long long hi) {
    const unsigned long long w = (unsigned long long)hi - (unsigned long long)lo + 1;
    return (long long)((unsigned long long)lo + (w == 0 ? rng() : rng() % w));
}
[[maybe_unused]] static long long elapsed_ms(std::chrono::steady_clock::time_point a) {
    return (long long)std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count();
}

using ll = long long;

static FILE* file_with(const std::string& s) {
    FILE* fp = std::tmpfile();
    std::fwrite(s.data(), 1, s.size(), fp);
    std::rewind(fp);
    return fp;
}
static std::string file_text(FILE* fp) {
    std::fflush(fp);
    std::rewind(fp);
    std::string s;
    char b[4096];
    size_t n;
    while ((n = std::fread(b, 1, sizeof b, fp)) > 0) s.append(b, n);
    return s;
}

static void random_test() {
    std::mt19937_64 rng(20261031);
    long long cases = 0;
    for (int it = 0; it < 300 * TEST_SCALE; ++it) {
        // 入力：整数・文字・文字列・実数を空白・改行・タブ混じりで
        std::ostringstream os;
        const int n = (int)rnd(rng, 0, 200);
        std::vector<int> kind(n);
        const char* sp[] = {" ", "\n", "  ", "\t", "\r\n"};
        for (int i = 0; i < n; ++i) {
            kind[i] = (int)rnd(rng, 0, 4);
            if (kind[i] == 0) os << (ll)rng();
            else if (kind[i] == 1) os << (char)('!' + rnd(rng, 0, 90));
            else if (kind[i] == 2) {
                const int L = (int)rnd(rng, 1, 8);
                for (int k = 0; k < L; ++k) os << (char)('a' + rnd(rng, 0, 25));
            } else if (kind[i] == 3) os << std::setprecision(17) << (double)(ll)rnd(rng, -1000000, 1000000) / 1024.0;
            else os << (unsigned long long)rng();
            os << sp[rnd(rng, 0, 4)];
        }
        const std::string text = os.str();
        FILE* fp = file_with(text);
        FastInput in(fp);
        std::istringstream is(text);
        for (int i = 0; i < n; ++i) {
            if (kind[i] == 0) {
                ll a, b;
                is >> a, in >> b;
                FI_CHECK(a == b);
            } else if (kind[i] == 1) {
                char a, b;
                is >> a, b = in.read<char>();
                FI_CHECK(a == b);
            } else if (kind[i] == 2) {
                std::string a, b;
                is >> a, in.read(b);
                FI_CHECK(a == b);
            } else if (kind[i] == 3) {
                double a, b;
                is >> a, b = in.read<double>();
                FI_CHECK(a == b);
            } else {
                unsigned long long a, b;
                is >> a, b = in.read<unsigned long long>();
                FI_CHECK(a == b);
            }
            ++cases;
        }
        FI_CHECK(in.eof());
        std::fclose(fp);
    }
    // 端の値・__int128
    {
        FILE* fp = file_with("-9223372036854775808 9223372036854775807 18446744073709551615 -170141183460469231731687303715884105728 +5 0 -0\n");
        FastInput in(fp);
        FI_CHECK(in.read<ll>() == std::numeric_limits<ll>::min());
        FI_CHECK(in.read<ll>() == std::numeric_limits<ll>::max());
        FI_CHECK(in.read<unsigned long long>() == std::numeric_limits<unsigned long long>::max());
        const __int128 m = in.read<__int128>();
        FI_CHECK(m < 0 && (unsigned __int128)m == (unsigned __int128)1 << 127);
        FI_CHECK(in.read<int>() == 5 && in.read<int>() == 0 && in.read<int>() == 0);
        FI_CHECK(in.eof());
        std::fclose(fp);
        cases += 6;
    }
    // 出力
    for (int it = 0; it < 200 * TEST_SCALE; ++it) {
        FILE* fp = std::tmpfile();
        std::ostringstream exp;
        {
            FastOutput out(fp);
            out.set_precision(6);
            exp << std::fixed << std::setprecision(6);
            const int n = (int)rnd(rng, 0, 300);
            for (int i = 0; i < n; ++i) {
                const int k = (int)rnd(rng, 0, 5);
                if (k == 0) {
                    const ll v = (ll)rng();
                    out << v << ' ', exp << v << ' ';
                } else if (k == 1) {
                    const int v = (int)rnd(rng, -1000, 1000);
                    out.print(v, std::string("ab"), 'c'), exp << v << " ab c\n";
                } else if (k == 2) {
                    const double v = (double)rnd(rng, -100000, 100000) / 7.0;
                    out << v << '\n', exp << v << '\n';
                } else if (k == 3) {
                    std::vector<int> v(rnd(rng, 0, 5));
                    for (auto& x : v) x = (int)rnd(rng, -9, 9);
                    out.print(v), exp << [&] {
                        std::string s;
                        for (size_t q = 0; q < v.size(); ++q) s += (q ? " " : "") + std::to_string(v[q]);
                        return s;
                    }() << '\n';
                } else if (k == 4) {
                    const unsigned long long v = rng();
                    out.write(v), exp << v;
                } else {
                    out << "xyz" << '\n', exp << "xyz\n";
                }
            }
            out.flush();
        }
        FI_CHECK(file_text(fp) == exp.str());
        std::fclose(fp);
        ++cases;
    }
    {
        FILE* fp = std::tmpfile();
        {
            FastOutput out(fp);
            out << std::numeric_limits<ll>::min() << ' ' << (__int128)(-1) * ((__int128)1 << 100) << '\n';
        }  // デストラクタで書き出す
        FI_CHECK(file_text(fp) == "-9223372036854775808 -1267650600228229401496703205376\n");
        std::fclose(fp);
        ++cases;
    }
    std::cout << "[tool_fast_io] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[tool_fast_io] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    const int N = TEST_SCALE == 1 ? 1000000 : 10000000;
    FILE* fp = std::tmpfile();
    auto t0 = std::chrono::steady_clock::now();
    {
        FastOutput out(fp);
        for (int i = 0; i < N; ++i) out << (ll)(rng() % 2000000001) - 1000000000 << '\n';
    }
    const ll tw = elapsed_ms(t0);
    std::rewind(fp);
    t0 = std::chrono::steady_clock::now();
    ll s = 0;
    {
        FastInput in(fp);
        for (int i = 0; i < N; ++i) s += in.read<ll>();
    }
    const ll tr = elapsed_ms(t0);
    // 比べる：stringstream（cin と同じ読み方）
    const std::string text = file_text(fp);
    t0 = std::chrono::steady_clock::now();
    ll s2 = 0;
    {
        std::istringstream is(text);
        for (int i = 0; i < N / 10; ++i) {
            ll x;
            is >> x, s2 += x;
        }
    }
    const ll ts = elapsed_ms(t0) * 10;
    std::fclose(fp);
    std::cout << "[tool_fast_io] benchmark (整数 " << N << " 個): 書く " << tw << " ms / 読む " << tr << " ms（istream >> だと約 " << ts
              << " ms）  (check " << (s + s2) % 1000 << ")\n";
}

}  // namespace fastio_test

// =============================================================================
//  使い方デモ
// =============================================================================
namespace fastio_demo {

[[maybe_unused]] static std::string to_s(bool b) { return b ? "true" : "false"; }
[[maybe_unused]] static std::string to_s(const std::string& s) { return "\"" + s + "\""; }
[[maybe_unused]] static std::string to_s(const char* s) { return to_s(std::string(s)); }
[[maybe_unused]] static std::string to_s(char c) { return std::string("'") + c + "'"; }
[[maybe_unused]] static std::string to_s(double v) {
    std::ostringstream os;
    os << std::setprecision(12) << v;
    return os.str();
}
template <class T>
static std::string to_s(const T& v);
template <class A, class B>
static std::string to_s(const std::pair<A, B>& p);
template <class T>
static std::string to_s(const std::vector<T>& v);
template <class T>
static std::string to_s(const std::optional<T>& v);
template <class T>
static std::string to_s(const T& v) {  // << で出せるもの（整数など）
    std::ostringstream os;
    os << v;
    return os.str();
}
template <class A, class B>
static std::string to_s(const std::pair<A, B>& p) {
    return "(" + to_s(p.first) + ", " + to_s(p.second) + ")";
}
template <class T>
static std::string to_s(const std::vector<T>& v) {
    std::string s = "[";
    for (size_t i = 0; i < v.size(); ++i) s += (i ? ", " : "") + to_s(v[i]);
    return s + "]";
}
template <class T>
static std::string to_s(const std::optional<T>& v) {
    return v ? to_s(*v) : std::string("nullopt");
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
[[maybe_unused]] static void note(const std::string& s) { std::cout << "    " << s << "\n"; }
[[maybe_unused]] static void code(const std::string& s) { std::cout << "  " << s << "\n"; }
[[maybe_unused]] static void stmt(const std::string& s) { code(s + ";"); }
template <class V>
[[maybe_unused]] static void show(const char* expr, const V& v, const std::string& comment) {
    std::string line = "  " + pad(expr, 38) + " -> " + to_s(v);
    if (!comment.empty()) line = pad(line, 66) + "  // " + comment;
    std::cout << line << "\n";
}
#define SHOW(expr, comment) fastio_demo::show(#expr, (expr), comment)
#define RUN(...)                  \
    fastio_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("高速入出力（tool_fast_io） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  fast_in: read<T>() / read(x, y, …) / >> / read_vec<T>(n) / eof()\n"
              << "  fast_out: write(x) / print(a, b, …) / << / set_precision(k) / flush()      FastInput(fp) / FastOutput(fp)\n";

    api("【1】 fast_in.read<T>() / fast_in.read(x, y, …) / fast_in >> x", "T / void / FastInput&", "標準入力から読む（cin の代わり。cin と混ぜない）",
        "O(文字数)");
    code("int N = fast_in.read<int>();");
    code("long long a, b; std::string s; fast_in.read(a, b, s);");
    code("auto A = fast_in.read_vec<long long>(N);");
    note("（デモは標準入力を使わないので、ファイルから読む形で見せる）");
    {
        FILE* fp = std::tmpfile();
        const char* text = "3\n10 -20 30\nhello 2.5 X\n";
        std::fwrite(text, 1, std::strlen(text), fp);
        std::rewind(fp);
        FastInput in(fp);
        RUN(const int N = in.read<int>());
        SHOW(N, "");
        SHOW(in.read_vec<long long>(N), "");
        RUN(std::string s; double d; char c);
        RUN(in.read(s, d, c));
        SHOW(s, "");
        SHOW(d, "");
        SHOW(c, "");
        SHOW(in.eof(), "もう読むものが無い");
        std::fclose(fp);
    }

    api("【2】 fast_out.print(a, b, …) / fast_out << x / fast_out.write(x)", "void / FastOutput& / void",
        "空白区切りで書いて改行 / cout の代わり（'\\n' も書く）/ 区切りなしで 1 つ", "O(文字数)");
    code("fast_out.print(N, ans);              // \"N ans\\n\"");
    code("fast_out.print(vec);                 // 中身を空白区切りで 1 行");
    code("fast_out << x << ' ' << y << '\\n';");
    code("fast_out.set_precision(10); fast_out << 3.14159265358979 << '\\n';");
    note("書いた分はたまっていて、プログラムの終わりに自動で書き出す（途中で cout を使うなら先に fast_out.flush()）");
}

}  // namespace fastio_demo

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
    fastio_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    fastio_test::benchmark();    // 速度計測
    fastio_demo::demo();         // 使い方デモ
    return 0;
}
