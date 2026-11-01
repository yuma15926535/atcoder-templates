// =============================================================================
//  インタラクティブ問題の補助（質問して返事を読む・flush・質問回数・手元のジャッジにつないで試す）   --- AtCoder 用テンプレート
// -----------------------------------------------------------------------------
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）
//  文字コード     : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  インタラクティブ問題では、質問を出すたびに flush しないと返事が来ない（cin.tie(nullptr) にしていると自動では流れない）。
//  Interactive は「1 行出して flush → 返事を読む」をまとめ、質問回数を数える。
//  同じ solve を、本番は cin / cout で、手元では自分で書いたジャッジ（関数）につないで何度でも試せる
//    例: ABC269 E（Last Rook）… io.ask("?", A, B, C, D) で数を聞き、二分探索 → io.answer("!", X, Y)
//        ABC244 C（Yamanote Line Game）… io.say(x) と io.read() を交互に
//
//  機能一覧                                                                                  計算量
//    Interactive io;                     本番：cin / cout でやり取りする
//    Interactive io(judge, first);       手元：judge は string(const string& 質問の 1 行) の関数で、返事の 1 行を返す。
//                                        first はジャッジが最初に送ってくる入力（N など。無ければ ""）
//    io.say(a, b, …)                     a b … を空白区切りの 1 行で出して flush（vector は中身を並べる）     O(出す量)
//    io.read<T>() / io.read_vec<T>(n)    返事・入力を 1 つ / n 個読む（T の既定は long long。string も可）
//    io.read_line()                      1 行まるごと読む
//    io.ask<R>(a, b, …)                  say(a, b, …) して read<R>() を返す。質問回数 +1
//    io.answer(a, b, …)                  答えを 1 行出す（質問回数に入れない）。answered() が true に
//    io.query_count()                    今までの質問回数
//    io.set_limit(q)                     質問回数の上限。超えて ask すると止める（手元で回数オーバーに気づける）
//    io.set_exit_on(v)                   読んだ値が v（例: −1 = ジャッジが不正と判断）ならすぐ exit(0)
//    io.set_log(os)                      やり取りを os（&std::cerr など）に表示する（nullptr で消す）
//    io.answered() / io.last_answer()    answer したか / 最後に answer した 1 行
//    interactive_tokens(line)            ジャッジ用：1 行を空白で区切った文字列の一覧
//    interactive_numbers(line)           ジャッジ用：1 行を空白で区切った long long の一覧（先頭の "?" などは飛ばす）
//
//  使い方
//    void solve(Interactive& io) { int N = io.read<int>(); … io.ask("?", l, r) … io.answer("!", x); }
//    本番：int main() { Interactive io; solve(io); }
//    手元：隠した答えを持つ judge を書き、Interactive io(judge, "N の行"); solve(io); のあと答えが合っているかを調べる
//  注意
//    say・ask・answer は毎回 flush する（endl と同じ）。Main_base の BUFF（cin.tie(nullptr)）のままでよい
//    手元のジャッジで、返事が足りないのに read すると止める（solve が余計に読もうとしている）
//  速さの目安: 手元のジャッジで ask を 10^5 回 約 0.05 秒
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    interactive_test::random_test()  … 数当て・ABC269 E の形のジャッジで何度も試す、本番の入出力の形も確認
//    interactive_test::benchmark()    … 速度計測
//    interactive_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
// =============================================================================
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <functional>
#include <iomanip>
#include <iostream>
#include <optional>
#include <random>
#include <sstream>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

// ======================= ここからテンプレート本体（コピー範囲） =======================
namespace interactive_internal {
[[noreturn]] inline void fail(const char* func, const char* msg) {
    std::cerr << func << ": " << msg << "\n";
    std::abort();
}
template <class T>
struct is_vec : std::false_type {};
template <class T, class A>
struct is_vec<std::vector<T, A>> : std::true_type {};
}  // namespace interactive_internal

struct Interactive {
    using Judge = std::function<std::string(const std::string&)>;
    Interactive() {}
    explicit Interactive(Judge judge, const std::string& first = "") : judge_(std::move(judge)), local_(true) { buf_.str(first + "\n"); }

    template <class... Args>
    void say(const Args&... args) {
        std::ostringstream os;
        bool first = true;
        (put(os, args, first), ...);
        (void)first;
        send(os.str());
    }
    template <class T = long long>
    T read() {
        T v{};
        if (local_) {
            if (!(buf_ >> v)) interactive_internal::fail("Interactive::read", "ジャッジからの返事が足りない（余計に読もうとしている）");
        } else {
            if (!(std::cin >> v)) interactive_internal::fail("Interactive::read", "入力が足りない");
        }
        if (log_) *log_ << "< " << v << "\n";
        if constexpr (std::is_arithmetic_v<T>)
            if (exit_on_ && (long long)v == *exit_on_) std::exit(0);
        return v;
    }
    template <class T = long long>
    std::vector<T> read_vec(int n) {
        std::vector<T> v(n);
        for (auto& x : v) x = read<T>();
        return v;
    }
    std::string read_line() {
        std::string s;
        std::istream& in = local_ ? static_cast<std::istream&>(buf_) : std::cin;
        in >> std::ws;
        if (!std::getline(in, s)) interactive_internal::fail("Interactive::read_line", "返事・入力が足りない");
        if (log_) *log_ << "< " << s << "\n";
        return s;
    }
    template <class R = long long, class... Args>
    R ask(const Args&... args) {
        if (limit_ >= 0 && queries_ >= limit_) interactive_internal::fail("Interactive::ask", "質問回数の上限を超えた");
        ++queries_;
        say(args...);
        return read<R>();
    }
    template <class... Args>
    void answer(const Args&... args) {
        std::ostringstream os;
        bool first = true;
        (put(os, args, first), ...);
        (void)first;
        last_answer_ = os.str();
        answered_ = true;
        send(last_answer_);
    }
    int query_count() const { return queries_; }
    void set_limit(int q) { limit_ = q; }
    void set_exit_on(long long v) { exit_on_ = v; }
    void set_log(std::ostream* os) { log_ = os; }
    bool answered() const { return answered_; }
    const std::string& last_answer() const { return last_answer_; }

   private:
    Judge judge_;
    bool local_ = false;
    std::stringstream buf_;  // 手元：ジャッジからの返事をためる
    int queries_ = 0, limit_ = -1;
    std::optional<long long> exit_on_;
    std::ostream* log_ = nullptr;
    bool answered_ = false;
    std::string last_answer_;

    template <class T>
    static void put(std::ostringstream& os, const T& x, bool& first) {
        if constexpr (interactive_internal::is_vec<T>::value) {
            for (const auto& e : x) put(os, e, first);
        } else {
            if (!first) os << ' ';
            os << x;
            first = false;
        }
    }
    void send(const std::string& line) {
        if (log_) *log_ << "> " << line << "\n";
        if (local_) {
            const std::string reply = judge_(line);
            // 読み残しの後ろに返事をつなぐ
            std::string rest;
            std::getline(buf_, rest, '\0');
            buf_.clear();
            buf_.str(rest + "\n" + reply + "\n");
        } else {
            std::cout << line << std::endl;  // flush する
        }
    }
};

// ---- ジャッジを書くとき用 ----
inline std::vector<std::string> interactive_tokens(const std::string& line) {
    std::istringstream is(line);
    std::vector<std::string> v;
    std::string t;
    while (is >> t) v.push_back(t);
    return v;
}
inline std::vector<long long> interactive_numbers(const std::string& line) {  // 数でないもの（"?" "!" など）は飛ばす
    std::vector<long long> v;
    for (const auto& t : interactive_tokens(line)) {
        char* end = nullptr;
        const long long x = std::strtoll(t.c_str(), &end, 10);
        if (end && *end == '\0' && !t.empty()) v.push_back(x);
    }
    return v;
}
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト・速度計測
// =============================================================================
namespace interactive_test {

#if defined(_DEBUG) || defined(_GLIBCXX_DEBUG) || (defined(__GNUC__) && !defined(__OPTIMIZE__))
[[maybe_unused]] static const int TEST_SCALE = 1;
#else
[[maybe_unused]] static const int TEST_SCALE = 8;
#endif

#define IA_CHECK(...)                                                                          \
    do {                                                                                       \
        if (!(__VA_ARGS__)) {                                                                  \
            std::cerr << "[tool_interactive] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
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

// 数当て：1〜N の隠した数 X。"? y" に "<"（X < y）/ ">="（X ≥ y）と答える。"! x" で答える
static ll guess_solve(Interactive& io) {
    const ll N = io.read();
    ll lo = 1, hi = N + 1;  // X ∈ [lo, hi)
    while (hi - lo > 1) {
        const ll mid = (lo + hi) / 2;
        if (io.ask<std::string>("?", mid) == "<") hi = mid;
        else lo = mid;
    }
    io.answer("!", lo);
    return lo;
}

// ABC269 E（Last Rook）の形：N × N に N − 1 個のルーク（どの 2 つも同じ行・列にない）。空いている行 X・列 Y を当てる
// "? A B C D" に [A, B] 行 × [C, D] 列のルークの数を答える（1-indexed・閉区間）
static std::pair<ll, ll> rook_solve(Interactive& io) {
    const ll N = io.read();
    ll lo = 1, hi = N;  // 空いている行は [lo, hi]
    while (lo < hi) {
        const ll mid = (lo + hi) / 2;
        if (io.ask("?", lo, mid, 1, N) < mid - lo + 1) hi = mid;
        else lo = mid + 1;
    }
    const ll X = lo;
    lo = 1, hi = N;
    while (lo < hi) {
        const ll mid = (lo + hi) / 2;
        if (io.ask("?", 1, N, lo, mid) < mid - lo + 1) hi = mid;
        else lo = mid + 1;
    }
    io.answer("!", X, lo);
    return {X, lo};
}

static void random_test() {
    std::mt19937_64 rng(20261031);
    long long cases = 0;
    // 数当て
    for (int it = 0; it < 2000 * TEST_SCALE; ++it) {
        const ll N = rnd(rng, 1, rnd(rng, 0, 1) ? 10 : 1000000000), X = rnd(rng, 1, N);
        ll lim = 0;
        while ((1LL << lim) < N) ++lim;
        bool ok = false;
        int questions = 0;
        auto judge = [&](const std::string& line) -> std::string {
            const auto t = interactive_tokens(line);
            IA_CHECK(!t.empty() && (t[0] == "?" || t[0] == "!"));
            if (t[0] == "!") {
                ok = std::stoll(t[1]) == X;
                return "";
            }
            ++questions;
            return X < std::stoll(t[1]) ? "<" : ">=";
        };
        Interactive io(judge, std::to_string(N));
        io.set_limit((int)lim);
        const ll got = guess_solve(io);
        IA_CHECK(ok && got == X && io.answered() && io.last_answer() == "! " + std::to_string(X));
        IA_CHECK(io.query_count() == questions && questions <= lim);
        ++cases;
    }
    // ABC269 E の形
    for (int it = 0; it < 500 * TEST_SCALE; ++it) {
        const int N = (int)rnd(rng, 2, rnd(rng, 0, 1) ? 8 : 1000);
        std::vector<int> col(N + 1);  // 行 r のルークの列（空いている行は 0）
        std::vector<int> cs(N);
        for (int i = 0; i < N; ++i) cs[i] = i + 1;
        std::shuffle(cs.begin(), cs.end(), rng);
        const int X = (int)rnd(rng, 1, N), Y = cs[N - 1];
        int k = 0;
        for (int r = 1; r <= N; ++r) col[r] = r == X ? 0 : cs[k++];
        std::pair<ll, ll> said{-1, -1};
        auto judge = [&](const std::string& line) -> std::string {
            const auto v = interactive_numbers(line);
            if (line[0] == '!') {
                said = {v[0], v[1]};
                return "";
            }
            IA_CHECK(v.size() == 4);
            int c = 0;
            for (int r = (int)v[0]; r <= v[1]; ++r) c += col[r] != 0 && v[2] <= col[r] && col[r] <= v[3];
            return std::to_string(c);
        };
        Interactive io(judge, std::to_string(N));
        io.set_limit(20);
        const auto got = rook_solve(io);
        IA_CHECK(got == std::make_pair((ll)X, (ll)Y) && said == got);
        ++cases;
    }
    // say：vector・文字・文字列が空白区切りの 1 行になる / read_vec / read_line / 読み残しの後ろに返事がつながる
    {
        std::vector<std::string> lines;
        auto judge = [&](const std::string& line) -> std::string {
            lines.push_back(line);
            return lines.size() == 1 ? "10 20 30" : "hello world";
        };
        Interactive io(judge, "3\nabc def");
        IA_CHECK(io.read<int>() == 3);
        io.say("?", std::vector<int>{1, 2, 3}, 'x', std::string("yz"), 4.5);
        IA_CHECK(lines.back() == "? 1 2 3 x yz 4.5");
        IA_CHECK(io.read_line() == "abc def");  // 最初の入力の残りが先
        IA_CHECK(io.read_vec<int>(3) == std::vector<int>({10, 20, 30}));
        io.say();
        IA_CHECK(lines.back().empty());
        IA_CHECK(io.read_line() == "hello world");
        IA_CHECK(interactive_numbers("? 1 -2 x 30") == std::vector<long long>({1, -2, 30}));
        cases += 6;
    }
    // 本番の形（cin / cout）：cin に返事を流し込み、cout に出たものを確かめる
    {
        std::istringstream fake_in("1000\n>=\n<\n>=\n<\n<\n>=\n>=\n>=\n<\n>=\n");
        std::ostringstream fake_out;
        auto* oi = std::cin.rdbuf(fake_in.rdbuf());
        auto* oo = std::cout.rdbuf(fake_out.rdbuf());
        Interactive io;
        const ll got = guess_solve(io);
        std::cin.rdbuf(oi), std::cout.rdbuf(oo);
        IA_CHECK(io.query_count() == 10 && got >= 1 && got <= 1000);
        const auto out = fake_out.str();
        IA_CHECK(out.substr(0, 6) == "? 501\n" && out.find("! " + std::to_string(got) + "\n") != std::string::npos);
        cases += 2;
    }
    std::cout << "[tool_interactive] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    if (TEST_SCALE == 1) std::cout << "[tool_interactive] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    const int Q = TEST_SCALE == 1 ? 10000 : 100000;
    ll s = 0;
    auto judge = [&](const std::string& line) { return std::to_string(interactive_numbers(line)[0] * 2); };
    const auto t0 = std::chrono::steady_clock::now();
    Interactive io(judge);
    for (int i = 0; i < Q; ++i) s += io.ask("?", i);
    std::cout << "[tool_interactive] benchmark: 手元のジャッジで ask を " << Q << " 回 " << elapsed_ms(t0) << " ms  (check " << s % 1000 << ")\n";
}

}  // namespace interactive_test

// =============================================================================
//  使い方デモ
// =============================================================================
namespace interactive_demo {

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
#define SHOW(expr, comment) interactive_demo::show(#expr, (expr), comment)
#define RUN(...)                  \
    interactive_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("インタラクティブ問題の補助（tool_interactive） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  Interactive io; / Interactive io(judge, first);\n"
              << "  say / read<T> / read_vec<T> / read_line / ask<R> / answer / query_count / set_limit / set_exit_on / set_log\n"
              << "  answered / last_answer      ジャッジ用: interactive_tokens / interactive_numbers\n";

    api("【1】 Interactive io(judge, first)：手元のジャッジにつなぐ", "—",
        "judge(質問の 1 行) が返事の 1 行を返す関数。first はジャッジが最初に送る入力。本番は Interactive io; だけ", "—");
    code("// 隠した数 X = 37（1〜100）。\"? y\" に X < y なら \"<\"、そうでなければ \">=\" と答えるジャッジ");
    RUN(const long long X = 37);
    RUN(auto judge = [&](const std::string& line) -> std::string { auto v = interactive_numbers(line); if (line[0] == '!') return ""; return X < v[0] ? "<" : ">="; });
    RUN(Interactive io(judge, "100"));

    api("【2】 read<T>() / ask<R>(…) / answer(…)", "T / R / void", "入力を読む / 1 行出して flush し返事を読む（質問回数 +1）/ 答えを出す", "—");
    RUN(io.set_log(&std::cout));  // やり取りを表示（本番では使わない）
    RUN(const long long N = io.read());
    RUN(long long lo = 1, hi = N + 1);
    RUN(while (hi - lo > 1) { long long mid = (lo + hi) / 2; if (io.ask<std::string>("?", mid) == "<") hi = mid; else lo = mid; });
    RUN(io.answer("!", lo));
    RUN(io.set_log(nullptr));
    SHOW(io.query_count(), "質問回数（⌈log2 100⌉ = 7 以下）");
    SHOW(io.answered(), "");
    SHOW(io.last_answer(), "");

    api("【3】 say(…) / read_vec<T>(n) / read_line()", "void / vector<T> / string", "1 行出す（vector は中身を並べる）/ n 個読む / 1 行読む", "—");
    RUN(auto echo = [](const std::string& line) { return "got " + line; });  // 質問をそのまま返すジャッジ
    RUN(Interactive e(echo, "3 1 4 1"));
    SHOW(e.read_vec<int>(4), "最初の入力");
    RUN(e.say("?", std::vector<int>{5, 9, 2}));
    SHOW(e.read_line(), "ジャッジの返事 1 行");

    api("【4】 set_limit(q) / set_exit_on(v)", "void", "質問回数の上限（超えて ask すると止める）/ 読んだ値が v ならすぐ終わる（−1 が来たら終了する問題用）", "—");
    code("io.set_limit(20);      // ABC269 E は 20 回まで");
    code("io.set_exit_on(-1);    // 不正な質問で −1 が返ってきたら exit(0)");

    api("【5】 本番での使い方", "—", "solve(Interactive&) を書いておけば、本番と手元のジャッジで同じものが動く", "—");
    code("void solve(Interactive& io) { … }");
    code("int main() { BUFF; Interactive io; solve(io); }               // 本番（cin / cout。毎回 flush する）");
    code("// 手元: for (試す回数) { 隠した答えを決めて judge を作る; Interactive io(judge, first); solve(io); 答えを確かめる; }");
}

}  // namespace interactive_demo

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
    interactive_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    interactive_test::benchmark();    // 速度計測
    interactive_demo::demo();         // 使い方デモ
    return 0;
}
