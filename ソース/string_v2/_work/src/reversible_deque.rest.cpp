// =============================================================================
//  ランダムテスト
// =============================================================================
namespace rdq_test {

{COMMON_TEST}

template <class T, class Gen>
static void test_deque(std::mt19937_64& rng, Gen gen, int& cases) {
    for (int it = 0; it < 400 * TEST_SCALE; ++it) {
        std::deque<T> ref;
        ReversibleDeque<T> d;
        if (rng() % 2) {
            std::vector<T> init(rng() % 20);
            for (auto& x : init) x = gen();
            d = ReversibleDeque<T>(init);
            ref.assign(init.begin(), init.end());
        }
        for (int q = 0; q < 300; ++q) {
            const int op = (int)(rng() % 12);
            if (op == 0) {
                const T x = gen();
                d.push_back(x), ref.push_back(x);
            } else if (op == 1) {
                const T x = gen();
                d.push_front(x), ref.push_front(x);
            } else if (op == 2 || op == 3) {
                std::vector<T> s(rng() % 40);
                for (auto& x : s) x = gen();
                if (op == 2) d.push_back(s), ref.insert(ref.end(), s.begin(), s.end());
                else d.push_front(s), ref.insert(ref.begin(), s.begin(), s.end());
            } else if (op == 4 && !ref.empty()) {
                RDQ_CHECK(d.pop_back() == ref.back());
                ref.pop_back();
            } else if (op == 5 && !ref.empty()) {
                RDQ_CHECK(d.pop_front() == ref.front());
                ref.pop_front();
            } else if (op == 6 || op == 7) {
                d.reverse(), std::reverse(ref.begin(), ref.end());
            } else if (op == 8 && !ref.empty()) {
                const int i = (int)(rng() % ref.size());
                const T x = gen();
                d[i] = x, ref[i] = x;
            } else if (op == 9 && rng() % 20 == 0) {
                d.clear(), ref.clear();
            } else if (op == 10 && rng() % 5 == 0) {
                d.reserve((int)(rng() % 200));  // 置き場を先に取っても中身は変わらない
            }
            RDQ_CHECK(d.size() == (int)ref.size() && d.empty() == ref.empty());
            if (!ref.empty()) {
                RDQ_CHECK(d.front() == ref.front() && d.back() == ref.back());
                const int i = (int)(rng() % ref.size());
                RDQ_CHECK(d[i] == ref[i]);
            }
            if (q % 30 == 0) {
                RDQ_CHECK(d.to_vector() == std::vector<T>(ref.begin(), ref.end()));
                const int pos = (int)(rng() % (ref.size() + 1)), len = (int)(rng() % (ref.size() - pos + 1));
                const auto sub = d.substr(pos, len);
                RDQ_CHECK((int)sub.size() == len);
                for (int k = 0; k < len; ++k) RDQ_CHECK(sub[k] == ref[pos + k]);
            }
        }
        ++cases;
    }
}

template <class T, class Gen>
static void test_rotatable(std::mt19937_64& rng, Gen gen, int& cases) {
    for (int it = 0; it < 400 * TEST_SCALE; ++it) {
        std::vector<T> ref(rng() % 30);
        for (auto& x : ref) x = gen();
        RotatableArray<T> ra(ref);
        const int n = (int)ref.size();
        for (int q = 0; q < 200; ++q) {
            const int op = (int)(rng() % 5);
            if (op == 0) {
                const long long k = rng() % 3 == 0 ? rnd_ll(rng) : (long long)(rng() % 70) - 35;
                ra.rotate_left(k);
                if (n > 0) std::rotate(ref.begin(), ref.begin() + (int)(((k % n) + n) % n), ref.end());
            } else if (op == 1) {
                const long long k = rng() % 3 == 0 ? rnd_ll(rng) : (long long)(rng() % 70) - 35;
                ra.rotate_right(k);
                if (n > 0) std::rotate(ref.begin(), ref.begin() + (int)(((-(k % n)) % n + n) % n), ref.end());
            } else if (op == 2) {
                ra.reverse(), std::reverse(ref.begin(), ref.end());
            } else if (op == 3 && n > 0) {
                const int i = (int)(rng() % n);
                const T x = gen();
                ra.set(i, x), ref[i] = x;
            } else if (op == 4 && n > 0) {
                const int i = (int)(rng() % n);
                const T x = gen();
                ra[i] = x, ref[i] = x;
            }
            RDQ_CHECK(ra.size() == n);
            if (n > 0) {
                const int i = (int)(rng() % n);
                RDQ_CHECK(ra[i] == ref[i] && ra.get(i) == ref[i]);
            }
            if (q % 20 == 0) RDQ_CHECK(ra.to_vector() == ref);
        }
        ++cases;
    }
}

static void random_test() {
    std::mt19937_64 rng(20261012);
    int cases = 0;
    test_deque<char>(rng, [&] { return (char)('a' + rng() % 3); }, cases);
    test_deque<long long>(rng, [&] { return rnd_ll(rng); }, cases);
    test_rotatable<char>(rng, [&] { return (char)('a' + rng() % 3); }, cases);
    test_rotatable<int>(rng, [&] { return (int)(rng() % 100) - 50; }, cases);
    // 文字列として
    {
        ReversibleDeque<char> d("bc");
        d.push_front("xa"), d.push_back(std::string("de"));
        RDQ_CHECK(d.to_string() == "xabcde" && d.substr(1, 3) == "abc");
        d.reverse();
        RDQ_CHECK(d.to_string() == "edcbax" && d.substr(0, 2) == "ed");
        d.push_front("12");
        RDQ_CHECK(d.to_string() == "12edcbax");
        // 自分の要素を足す（置き場が広がる瞬間も）・自分自身をつなげる
        for (int len = 1; len <= 70; ++len) {
            std::string t;
            for (int i = 0; i < len; ++i) t += (char)('a' + i % 26);
            ReversibleDeque<char> e(t);
            e.push_back(e.front()), e.pop_front();
            RDQ_CHECK(e.to_string() == t.substr(1) + t[0]);
            ReversibleDeque<char> f(t);
            f.push_front(f.back()), f.pop_back();
            RDQ_CHECK(f.to_string() == t.back() + t.substr(0, t.size() - 1));
            ReversibleDeque<char> g(t);
            g.push_back(g), g.push_front(g);
            RDQ_CHECK(g.to_string() == t + t + t + t);
        }
        RotatableArray<char> ra(std::string("abcde"));
        ra.rotate_left(2);
        RDQ_CHECK(ra.to_string() == "cdeab");
        ra.reverse();
        RDQ_CHECK(ra.to_string() == "baedc");
        ra.rotate_right(1);
        RDQ_CHECK(ra.to_string() == "cbaed");
        ++cases;
    }
    std::cout << "[reversible_deque] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[reversible_deque] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    auto ms = [](std::chrono::steady_clock::time_point a) {
        return (long long)std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count();
    };
    const int Q = 10000000;
    auto t0 = std::chrono::steady_clock::now();
    ReversibleDeque<char> d;
    long long x = 0;
    for (int q = 0; q < Q; ++q) {
        const unsigned r = (unsigned)(rng() % 8);
        if (r == 0) d.reverse();
        else if (r <= 3) d.push_front((char)('a' + r));
        else if (r <= 6) d.push_back((char)('a' + r));
        else if (!d.empty()) x += d.pop_front();
    }
    x += (long long)d.to_string().size();
    const long long t_d = ms(t0);
    RotatableArray<int> ra(std::vector<int>(1000000, 1));
    t0 = std::chrono::steady_clock::now();
    for (int q = 0; q < Q; ++q) {
        const unsigned r = (unsigned)(rng() % 4);
        if (r == 0) ra.reverse();
        else if (r == 1) ra.rotate_left((long long)(rng() % 1000000));
        else x += ra[(int)(rng() % 1000000)];
    }
    const long long t_r = ms(t0);
    std::cout << "[reversible_deque] benchmark: ReversibleDeque 操作 10^7 回 " << t_d << " ms, RotatableArray（N = 10^6）操作 10^7 回 " << t_r << " ms  (check " << x % 7 << ")\n";
}

}  // namespace rdq_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace rdq_demo {

{COMMON_DEMO}

static void demo() {
    banner("前後に足せて反転できる列（reversible_deque） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  ReversibleDeque<T>: push_back / push_front（1 つ・列）/ pop_back / pop_front / front / back / [] / reverse / size / to_vector / to_string / substr / reserve\n"
              << "  RotatableArray<T>:  rotate_left / rotate_right / reverse / [] / get / set / size / to_vector / to_string\n"
              << "  反転・回転は O(1)。T = char で文字列（to_string・substr が string）\n";

    api("【1】 ReversibleDeque<T> d(s)", "push_back / push_front / pop_back / pop_front / front / back / d[i]",
        "前後への追加・取り出しと i 番目。列を足すとそのままの順でつながる（push_front(s) は s + 今の列）", "O(1)（列は O(|s|)）");
    {
        RUN(ReversibleDeque<char> d("abc"));
        RUN(d.push_front('x'));
        RUN(d.push_back("yz"));
        SHOW(d.to_string(), "");
        RUN(d.push_front("12"));
        SHOW(d.to_string(), "");
        SHOW(d.pop_back(), "");
        SHOW(d.front(), "");
        SHOW(d[3], "");

        api("【2】 d.reverse()", "なし", "全体を反転（印を変えるだけ。この後の push_front は反転後の先頭に足す）", "O(1)");
        RUN(d.reverse());
        SHOW(d.to_string(), "");
        RUN(d.push_front('<'));
        RUN(d.push_back('>'));
        SHOW(d.to_string(), "");
        SHOW(d.substr(1, 3), "1 番目から 3 文字");
        SHOW(d.reversed(), "");
    }
    {
        note("ABC158 D（String Formation）: 反転・先頭に足す・末尾に足す をそのまま。入力例 1:");
        ReversibleDeque<char> d("a");
        d.push_front('p'), d.reverse(), d.push_back('c'), d.reverse();
        note("  a → 前に p → 反転 → 後ろに c → 反転 = " + d.to_string());
    }
    {
        note("ABC237 D（LR insertion）: 後ろから見る。[N] から始め、s_i が L なら i - 1 を末尾に、R なら先頭に足す。入力例 1（LRRLR）:");
        const std::string S = "LRRLR";
        const int N = (int)S.size();
        ReversibleDeque<int> d;
        d.push_back(N);
        for (int i = N; i >= 1; --i) {
            if (S[i - 1] == 'L') d.push_back(i - 1);
            else d.push_front(i - 1);
        }
        note("  → " + join(d.to_vector()));
    }

    api("【3】 RotatableArray<T> ra(s)", "rotate_left(k) / rotate_right(k) / reverse() / ra[i] / get / set",
        "大きさが決まった配列の回転（新しい [i] = 元の [i + k]）と反転を O(1)。k は負・N 以上でもよい", "O(1)");
    {
        RUN(RotatableArray<char> ra(std::string("abcdef")));
        RUN(ra.rotate_left(2));
        SHOW(ra.to_string(), "");
        RUN(ra.reverse());
        SHOW(ra.to_string(), "");
        RUN(ra.rotate_right(1));
        SHOW(ra.to_string(), "");
        RUN(ra.set(0, 'Z'));
        SHOW(ra.to_string(), "");
        SHOW(ra[1], "");
    }
    {
        note("ABC199 C（IPFL）: 前半と後半の入れ替え = rotate_left(N)。入力例 2:");
        const int N = 2;
        RotatableArray<char> ra(std::string("FLIP"));
        const int qs[][3] = {{1, 1, 3}, {2, 0, 0}, {1, 1, 2}, {1, 2, 3}, {2, 0, 0}, {1, 1, 4}};
        for (const auto& q : qs) {
            if (q[0] == 2) {
                ra.rotate_left(N);
            } else {
                const char a = ra[q[1] - 1];
                ra[q[1] - 1] = ra[q[2] - 1], ra[q[2] - 1] = a;
            }
        }
        note("  → " + ra.to_string());
    }
}

}  // namespace rdq_demo

{SETUP_CONSOLE}

int main() {
    setup_console();
    rdq_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    rdq_test::benchmark();    // 速度計測
    rdq_demo::demo();         // 使い方デモ
    return 0;
}
