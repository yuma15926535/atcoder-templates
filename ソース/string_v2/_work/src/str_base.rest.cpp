// =============================================================================
//  ランダムテスト
// =============================================================================
namespace strb_test {

{COMMON_TEST}

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

{COMMON_DEMO}

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

{SETUP_CONSOLE}

int main() {
    setup_console();
    strb_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    strb_test::benchmark();    // 速度計測
    strb_demo::demo();         // 使い方デモ
    return 0;
}
