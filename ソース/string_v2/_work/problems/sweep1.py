# 最終照合の広い洗い出し（AtCoder_文字列問題一覧.txt の問題。作ったものを見ていないエージェントが解いたもの）
# agent_1 の担当問題（string_v2 最終照合）
# 使い方: python3 /home/claude/final/run.py probs.py /mnt/user-data/uploads/claude/Main_base.cpp /home/claude/final/out /home/claude/ac-library [問題名…]
MATH = '/mnt/user-data/uploads/claude/math_v2/'

# ---------- 答えが複数ある問題・小数の問題の判定 ----------
def chk_280c(inp, out, exp):
    S, T = inp.split()
    a = int(out.split()[0])
    return 1 <= a <= len(T) and T[:a - 1] + T[a:] == S

def chk_459d(inp, out, exp):
    toks = inp.split()
    t = int(toks[0]); cases = toks[1:1 + t]
    o = out.split(); e = exp.split(); oi = ei = 0
    for s in cases:
        if oi >= len(o) or o[oi] != e[ei]: return False
        if e[ei] == 'Yes':
            r = o[oi + 1]; oi += 2; ei += 2
            if sorted(r) != sorted(s) or any(r[i] == r[i + 1] for i in range(len(r) - 1)): return False
        else:
            oi += 1; ei += 1
    return oi == len(o)

def chk_float(inp, out, exp):
    return abs(float(out.split()[0]) - float(exp.split()[0])) <= 1e-9

# ======================= C〜F =======================
P('ABC280_C', ['00_str_base'], r'''
void solv(){
    string S, T; cin >> S >> T;
    cout << common_prefix_length(S, T) + 1 << '\n';  // 最初に食い違う位置（1 始まり）
}''', [("atcoder\natcorder\n", "5\n"), ("million\nmilllion\n", "5\n"), ("vvwvw\nvvvwvw\n", "3\n")], check=chk_280c)

P('ABC394_C', ['00_str_base'], r'''
void solv(){
    string S; cin >> S;
    // 左端の WA を AC に、を繰り返す ＝ 右からなめて WA を AC にする（W…WA → AC…C）
    for (int i = (int)S.size() - 2; i >= 0; --i)
        if (S[i] == 'W' && S[i + 1] == 'A') S[i] = 'A', S[i + 1] = 'C';
    cout << S << '\n';
}''', [("WACWA\n", "ACCAC\n"), ("WWA\n", "ACC\n"), ("WWWWW\n", "WWWWW\n")])

P('ABC299_C', ['00_str_base'], r'''
void solv(){
    int N; string S; cin >> N >> S;
    const auto rl = run_length(S);
    long long ans = -1;
    if (rl.size() >= 2)  // o と - の両方がある → o の塊はどれも - と隣り合う
        for (const auto& b : rl) if (b.value == 'o') ans = max(ans, b.len);
    cout << ans << '\n';
}''', [("10\no-oooo---o\n", "4\n"), ("1\n-\n", "-1\n"), ("30\n-o-o-oooo-oo-o-ooooooo--oooo-o\n", "7\n")])

P('ABC380_C', ['00_str_base'], r'''
void solv(){
    int N, K; string S; cin >> N >> K >> S;
    auto rl = run_length(S);
    int cnt = 0;
    for (int i = 0; i < (int)rl.size(); ++i)
        if (rl[i].value == '1' && ++cnt == K) { swap(rl[i - 1], rl[i]); break; }  // 直前の 0 の塊と入れ替える
    cout << run_length_decode(rl) << '\n';
}''', [("15 3\n010011100011001\n", "010011111000001\n"), ("10 2\n1011111111\n", "1111111110\n")])

P('ABC433_C', ['00_str_base'], r'''
void solv(){
    string S; cin >> S;
    const auto rl = run_length(S);
    long long ans = 0;
    for (size_t i = 0; i + 1 < rl.size(); ++i)
        if (rl[i].value + 1 == rl[i + 1].value) ans += min(rl[i].len, rl[i + 1].len);
    cout << ans << '\n';
}''', [("1122\n", "2\n"), ("7788788\n", "3\n"), ("2025\n", "0\n"), ("1112222334445556555\n", "11\n")])

P('ABC372_C', ['02_pattern_match'], r'''
void solv(){
    int N, Q; string S; cin >> N >> Q >> S;
    long long cnt = count_occurrences(S, "ABC");
    auto at = [&](int i) { return i >= 0 && i + 3 <= N && S.compare(i, 3, "ABC") == 0; };
    while (Q--) {
        int X; char C; cin >> X >> C; --X;
        for (int i = X - 2; i <= X; ++i) cnt -= at(i);
        S[X] = C;
        for (int i = X - 2; i <= X; ++i) cnt += at(i);
        cout << cnt << '\n';
    }
}''', [("7 4\nABCDABC\n4 B\n3 A\n5 C\n4 G\n", "2\n1\n1\n0\n"),
       ("3 3\nABC\n1 A\n2 B\n3 C\n", "1\n1\n1\n"),
       ("15 10\nBBCCBCACCBACACA\n9 C\n11 B\n5 B\n11 B\n4 A\n8 C\n8 B\n5 B\n7 B\n14 B\n", "0\n0\n0\n0\n1\n1\n2\n2\n1\n1\n")])

P('ABC264_D', ['00_str_base', MATH + '04_permutation'], r'''
void solv(){
    string S; cin >> S;
    // "atcoder" の順の番号の列にして、転倒数 = 隣の入れ替えの最小回数
    cout << inversion_count(to_ids(S, Alphabet::of("atcoder"))) << '\n';
}''', [("catredo\n", "8\n"), ("atcoder\n", "0\n"), ("redocta\n", "21\n")])

P('ABC345_C', ['00_str_base'], r'''
void solv(){
    string S; cin >> S;
    const auto c = count_chars(S);
    const long long n = S.size();
    long long ans = n * (n - 1) / 2;
    bool dup = false;
    for (long long x : c) { ans -= x * (x - 1) / 2; if (x >= 2) dup = true; }  // auto x だと int で x*(x-1) が溢れる
    cout << ans + dup << '\n';
}''', [("abc\n", "3\n"), ("aaaaa\n", "1\n")])

P('ABC456_D', [MATH + '00_modint'], r'''
using mint = ModInt998;
void solv(){
    string S; cin >> S;
    array<mint, 3> dp{};  // dp[c] = 最後が c の部分列（隣が同じでない）の数
    for (char ch : S) {
        const int x = ch - 'a';
        mint add = 1;
        for (int c = 0; c < 3; ++c) if (c != x) add += dp[c];
        dp[x] += add;
    }
    cout << dp[0] + dp[1] + dp[2] << '\n';
}''', [("abbc\n", "11\n"),
       ("cabcabcbcaccacbcbcaabacbacaabccacbccbcacbacbacabcacabcaccaaaaabababcbabacaccabbcacbcbcbcababcbcbabca\n", "378217423\n")])

P('ABC459_D', ['00_str_base'], r'''
void solv(){
    int T; cin >> T;
    while (T--) {
        string S; cin >> S;
        const int n = S.size();
        const auto c = count_chars(S);
        if (*max_element(nall(c)) > (n + 1) / 2) { NO; continue; }
        vector<int> ord(26);
        iota(nall(ord), 0);
        sort(nall(ord), [&](int a, int b) { return c[a] > c[b]; });
        string g;  // 多い文字から順に並べたもの
        for (int k : ord) g += string(c[k], 'a' + k);
        string t(n, ' ');
        int p = 0;
        for (int i = 0; i < n; i += 2) t[i] = g[p++];  // 偶数番目 → 奇数番目の順に詰める
        for (int i = 1; i < n; i += 2) t[i] = g[p++];
        YES;
        cout << t << '\n';
    }
}''', [("3\naiiw\ndoodoo\naabbababcacababaaba\n", "Yes\niwai\nNo\nYes\nababacabababacababa\n")], check=chk_459d)

P('ABC287_D', ['00_str_base'], r'''
void solv(){
    string S, T; cin >> S >> T;
    const int n = S.size(), m = T.size();
    auto ok = [](char a, char b) { return a == '?' || b == '?' || a == b; };
    int P = 0, Q = 0;  // ? 込みで一致する接頭辞・接尾辞の長さ
    while (P < m && ok(S[P], T[P])) ++P;
    while (Q < m && ok(S[n - 1 - Q], T[m - 1 - Q])) ++Q;
    for (int x = 0; x <= m; ++x) { if (x <= P && m - x <= Q) YES; else NO; }
}''', [("a?c\nb?\n", "Yes\nNo\nNo\n"), ("atcoder\n?????\n", "Yes\nYes\nYes\nYes\nYes\nYes\n"),
       ("beginner\ncontest\n", "No\nNo\nNo\nNo\nNo\nNo\nNo\nNo\n")])

P('ABC344_D', ['02_pattern_match'], r'''
void solv(){
    string T; int N; cin >> T >> N;
    const int m = T.size();
    vector<int> dp(m + 1, INF);
    dp[0] = 0;
    rep(i, N) {
        int A; cin >> A;
        auto nd = dp;
        rep(a, A) {
            string s; cin >> s;
            for (int p : find_all(T, s)) if (dp[p] < INF) chmin(nd[p + s.size()], dp[p] + 1);
        }
        dp = nd;
    }
    cout << (dp[m] >= INF ? -1 : dp[m]) << '\n';
}''', [("abcde\n3\n3 ab abc abcd\n4 f c cd bcde\n2 e de\n", "2\n"),
       ("abcde\n3\n2 ab abc\n3 f c bcde\n1 e\n", "-1\n"),
       ("aaabbbbcccc\n6\n2 aa aaa\n2 dd ddd\n2 ab aabb\n4 bbaa bbbc bbb bbcc\n2 cc bcc\n3 ccc cccc ccccc\n", "4\n")])

P('ABC295_D', ['00_str_base'], r'''
void solv(){
    string S; cin >> S;
    vector<long long> cnt(1 << 10, 0);
    int mask = 0;
    long long ans = 0;
    cnt[0] = 1;
    for (char c : S) { mask ^= 1 << (c - '0'); ans += cnt[mask]++; }  // 各数字の個数の偶奇が同じ先頭どうしを組にする
    cout << ans << '\n';
}''', [("20230322\n", "4\n"), ("0112223333444445555556666666777777778888888889999999999\n", "185\n"),
       ("3141592653589793238462643383279502884197169399375105820974944\n", "9\n")])

P('ABC401_D', ['00_str_base'], r'''
void solv(){
    int N; long long K; string S; cin >> N >> K >> S;
    rep(i, N) if (S[i] == 'o') { if (i > 0) S[i - 1] = '.'; if (i + 1 < N) S[i + 1] = '.'; }
    const long long cnt = count(nall(S), 'o');
    if (cnt == K) {
        for (auto& c : S) if (c == '?') c = '.';
    } else {
        const auto rl = run_length(S);
        long long M = 0;  // ? の塊に置ける o の最大の数の和
        for (const auto& b : rl) if (b.value == '?') M += (b.len + 1) / 2;
        if (cnt + M == K)
            for (const auto& b : rl)
                if (b.value == '?' && b.len % 2 == 1) rep(j, b.len) S[b.start + j] = j % 2 ? '.' : 'o';
    }
    cout << S << '\n';
}''', [("4 2\no???\n", "o.??\n"), ("5 2\n?????\n", "?????\n"), ("7 3\n.o???o.\n", ".o.o.o.\n")])

P('ABC242_D', ['00_str_base'], r'''
void solv(){
    string S; int Q; cin >> S >> Q;
    const Alphabet A = Alphabet::of("ABC");
    while (Q--) {
        long long t, k; cin >> t >> k; --k;
        long long sh = 0;  // 親をたどる：k 文字目は k/2 文字目の子で、(k%2 + 1) ずれる
        while (t > 0 && k > 0) { sh += k % 2 + 1; k /= 2; --t; }
        sh += t;           // 先頭のまま残った t 回は 1 ずつずれる
        cout << A.ch((A.id(S[k]) + sh) % 3) << '\n';
    }
}''', [("ABC\n4\n0 1\n1 1\n1 3\n1 6\n", "A\nB\nC\nB\n"),
       ("CBBAACCCCC\n5\n57530144230160008 659279164847814847\n29622990657296329 861239705300265164\n509705228051901259 994708708957785197\n176678501072691541 655134104344481648\n827291290937314275 407121144297426665\n", "A\nA\nC\nA\nA\n")])

P('ABC242_E', ['00_str_base', MATH + '00_modint'], r'''
using mint = ModInt998;
void solv(){
    int T; cin >> T;
    while (T--) {
        int N; string S; cin >> N >> S;
        const int h = (N + 1) / 2;
        mint v = 0;  // 前半が S の前半より小さい回文の数 = 前半を 26 進数で読んだ値
        rep(i, h) v = v * 26 + (S[i] - 'A');
        const string P = S.substr(0, h) + reversed(S.substr(0, N / 2));  // 前半が同じ回文
        if (P <= S) v += 1;
        cout << v << '\n';
    }
}''', [("5\n3\nAXA\n6\nABCZAZ\n30\nQWERTYUIOPASDFGHJKLZXCVBNMQWER\n28\nJVIISNEOXHSNEAAENSHXOENSIIVJ\n31\nKVOHEEMSOZZASHENDIGOJRTJVMVSDWW\n",
        "24\n29\n212370247\n36523399\n231364016\n")])

# ======================= B =======================
P('ABC218_B', ['00_str_base'], r'''
void solv(){
    vector<int> P(26); cin >> P;   // vector<ll> だと from_ids に渡せない（コンパイルエラー）
    for (auto& x : P) --x;
    cout << from_ids(P) << '\n';
}''', [("1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26\n", "abcdefghijklmnopqrstuvwxyz\n"),
       ("2 1 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26\n", "bacdefghijklmnopqrstuvwxyz\n"),
       ("5 11 12 16 25 17 18 1 7 10 4 23 20 3 2 24 26 19 14 9 6 22 8 13 15 21\n", "eklpyqragjdwtcbxzsnifvhmou\n")])

P('ABC219_B', ['00_str_base'], r'''
void solv(){
    vector<string> S(3); cin >> S;
    string T; cin >> T;
    string ans;
    for (char c : T) ans += S[c - '1'];
    cout << ans << '\n';
}''', [("mari\nto\nzzo\n1321\n", "marizzotomari\n"), ("abra\ncad\nabra\n123\n", "abracadabra\n"), ("a\nb\nc\n1\n", "a\n")])

P('ABC382_B', ['00_str_base'], r'''
void solv(){
    int N, D; string S; cin >> N >> D >> S;
    for (int i = N - 1; i >= 0 && D > 0; --i) if (S[i] == '@') S[i] = '.', --D;
    cout << S << '\n';
}''', [("5 2\n.@@.@\n", ".@...\n"), ("3 3\n@@@\n", "...\n"), ("10 4\n@@@.@@.@@.\n", "@@@.......\n")])

P('ABC357_B', ['00_str_base'], r'''
void solv(){
    string S; cin >> S;
    const int up = count_if(nall(S), [](char c) { return isupper(c); });  // count_chars(S, Alphabet::upper()) は小文字で止まる
    cout << (up * 2 > (int)S.size() ? to_upper(S) : to_lower(S)) << '\n';
}''', [("AtCoder\n", "atcoder\n"), ("SunTORY\n", "SUNTORY\n"), ("a\n", "a\n")])

P('ABC476_B', ['02_wildcard_match'], r'''
void solv(){
    int N; string S, T; cin >> N >> S >> T;
    if (!wildcard_match(S, T, '*').empty()) YN
}''', [("8\nchokudai\n**o*u*ai\n", "Yes\n"), ("5\nsnuke\nsnake\n", "No\n"), ("5\nyiwiy\n*****\n", "Yes\n")])

P('ABC286_B', ['02_pattern_match'], r'''
void solv(){
    int N; string S; cin >> N >> S;
    cout << replace_all(S, "na", "nya") << '\n';
}''', [("4\nnaan\n", "nyaan\n"), ("4\nnear\n", "near\n"), ("8\nnational\n", "nyationyal\n")])

P('ABC414_B', ['00_str_base'], r'''
void solv(){
    int N; cin >> N;
    vector<RunBlock<char>> rl(N);
    long long tot = 0;
    for (auto& b : rl) { cin >> b.value >> b.len; tot = min(tot + b.len, 101LL); }  // 先に長さを確かめる（decode は 2×10^9 超で止まる）
    if (tot > 100) cout << "Too Long\n";
    else cout << run_length_decode(rl) << '\n';
}''', [("8\nm 1\ni 1\ns 2\ni 1\ns 2\ni 1\np 2\ni 1\n", "mississippi\n"),
       ("7\na 1000000000000000000\nt 1000000000000000000\nc 1000000000000000000\no 1000000000000000000\nd 1000000000000000000\ne 1000000000000000000\nr 1000000000000000000\n", "Too Long\n"),
       ("1\na 100\n", "a" * 100 + "\n"),
       ("6\ng 4\nj 1\nm 4\ne 4\nd 3\ni 4\n", "ggggjmmmmeeeedddiiii\n")])

P('ABC412_B', ['00_str_base'], r'''
void solv(){
    string S, T; cin >> S >> T;
    bool ok = true;
    for (size_t i = 1; i < S.size(); ++i)
        if (isupper(S[i]) && T.find(S[i - 1]) == string::npos) ok = false;  // Alphabet::of(T) は T に同じ文字があると止まる
    if (ok) YN
}''', [("AtCoder\nTotal\n", "Yes\n"), ("aBCdE\nabcdcba\n", "No\n"), ("abcde\nXYZ\n", "Yes\n")])

P('ABC242_B', ['00_str_base'], r'''
void solv(){
    string S; cin >> S;
    cout << sorted_chars(S) << '\n';
}''', [("aba\n", "aab\n"), ("zzzz\n", "zzzz\n")])

P('ABC297_B', ['00_str_base'], r'''
void solv(){
    string S; cin >> S;
    const CharPositions P(S, Alphabet::upper());
    const bool ok = (P.kth('B', 0) + P.kth('B', 1)) % 2 == 1 && P.kth('R', 0) < P.kth('K', 0) && P.kth('K', 0) < P.kth('R', 1);
    if (ok) YN
}''', [("RNBQKBNR\n", "Yes\n"), ("KRRBBNNQ\n", "No\n"), ("BRKRBQNN\n", "No\n")])

P('ABC447_B', ['00_str_base'], r'''
void solv(){
    string S; cin >> S;
    const auto c = count_chars(S);
    const int mx = *max_element(nall(c));
    string ans;
    for (char ch : S) if (c[ch - 'a'] != mx) ans += ch;
    cout << ans << '\n';
}''', [("mississippi\n", "mpp\n"), ("atcoder\n", "\n"), ("beginner\n", "bgir\n")])

P('ABC322_B', ['00_str_base'], r'''
void solv(){
    int N, M; string S, T; cin >> N >> M >> S >> T;
    const bool pre = common_prefix_length(S, T) == N, suf = common_suffix_length(S, T) == N;
    cout << (pre ? (suf ? 0 : 1) : (suf ? 2 : 3)) << '\n';
}''', [("3 7\nabc\nabcdefg\n", "1\n"), ("3 4\nabc\naabc\n", "2\n"), ("3 3\nabc\nxyz\n", "3\n"), ("3 3\naaa\naaa\n", "0\n")])

P('ABC352_B', ['06_subsequence'], r'''
void solv(){
    string S, T; cin >> S >> T;
    const NextTable nx(T);
    vector<int> ans;
    int p = 0;
    for (char c : S) { p = nx.next(p, c); ans.push_back(p + 1); ++p; }  // 貪欲に当てはめた位置
    cout << join(ans) << '\n';
}''', [("abc\naxbxyc\n", "1 3 6\n"), ("aaaa\nbbbbaaaa\n", "5 6 7 8\n"), ("atcoder\natcoder\n", "1 2 3 4 5 6 7\n")])

P('ABC221_B', ['00_str_base'], r'''
void solv(){
    string S, T; cin >> S >> T;
    bool ok = S == T;
    for (size_t i = 0; i + 1 < S.size(); ++i) { swap(S[i], S[i + 1]); ok |= S == T; swap(S[i], S[i + 1]); }
    if (ok) YN
}''', [("abc\nacb\n", "Yes\n"), ("aabb\nbbaa\n", "No\n"), ("abcde\nabcde\n", "Yes\n")])

P('ABC212_B', ['00_str_base'], r'''
void solv(){
    string X; cin >> X;
    const bool weak = run_length(X).size() == 1 || caesar_shift(X.substr(0, 3), 1, Alphabet::digits()) == X.substr(1);
    cout << (weak ? "Weak" : "Strong") << '\n';
}''', [("7777\n", "Weak\n"), ("0112\n", "Strong\n"), ("9012\n", "Weak\n")])

P('ABC441_B', ['00_str_base'], r'''
void solv(){
    int N, M; string S, T; cin >> N >> M >> S >> T;
    const Alphabet A = Alphabet::of(S), B = Alphabet::of(T);  // 文字の集合として（S・T は相異なる文字）
    int Q; cin >> Q;
    while (Q--) {
        string w; cin >> w;
        const bool inS = all_of(nall(w), [&](char c) { return A.contains(c); });
        const bool inT = all_of(nall(w), [&](char c) { return B.contains(c); });
        cout << (inS && !inT ? "Takahashi" : !inS && inT ? "Aoki" : "Unknown") << '\n';
    }
}''', [("6 5\nahikst\naikot\n5\nasahi\nokita\nkiai\nhash\nit\n", "Takahashi\nAoki\nUnknown\nTakahashi\nUnknown\n"),
       ("7 6\nahiknst\nahikos\n5\nkioki\nohiki\ntashi\nnishi\nkashi\n", "Aoki\nAoki\nTakahashi\nTakahashi\nUnknown\n"),
       ("13 11\ndefghiqsvwxyz\nacejmoqrtwx\n15\nqhsqzhd\njcareec\nwwqxqew\nwxqxwex\njxxrtwa\ntrtqjxe\nsqyggse\nxxqwxew\nxewwxxw\nwwqxwex\nxqqxqwq\nqxxexxe\nteqeroc\neeeqqee\nvxdevyy\n",
        "Takahashi\nAoki\nUnknown\nUnknown\nAoki\nAoki\nTakahashi\nUnknown\nUnknown\nUnknown\nUnknown\nUnknown\nAoki\nUnknown\nTakahashi\n")])

P('ABC347_B', ['03_suffix_array'], r'''
void solv(){
    string S; cin >> S;
    cout << SuffixArray(S).count_distinct_substrings() << '\n';
}''', [("yay\n", "5\n"), ("aababc\n", "17\n"), ("abracadabra\n", "54\n")])

P('ABC337_B', ['00_str_base'], r'''
void solv(){
    string S; cin >> S;
    if (is_sorted(nall(S))) YN
}''', [("AAABBBCCCCCCC\n", "Yes\n"), ("ACABABCBC\n", "No\n"), ("A\n", "Yes\n"), ("ABBBBBBBBBBBBBCCCCCC\n", "Yes\n")])

P('ABC418_B', ['00_str_base'], r'''
void solv(){
    string S; cin >> S;
    const int n = S.size();
    const CharPrefix cp(S);
    double best = 0;
    rep(l, n) rep(r, l + 3, n + 1)
        if (S[l] == 't' && S[r - 1] == 't') chmax(best, (double)(cp.count(l, r, 't') - 2) / (r - l - 2));
    cout << FIXED(17) << best << '\n';
}''', [("attitude\n", "0.50000000000000000\n"), ("ottottott\n", "0.66666666666666667\n"), ("coffeecup\n", "0.00000000000000000\n")], check=chk_float)

P('ABC346_B', ['00_str_base'], r'''
void solv(){
    int W, B; cin >> W >> B;
    const string S = repeat_string("wbwbwwbwbwbw", 20);
    const CharPrefix cp(S, Alphabet::of("wb"));
    bool ok = false;
    rep(l, 12) if (cp.count(l, l + W + B, 'w') == W && cp.count(l, l + W + B, 'b') == B) ok = true;
    if (ok) YN
}''', [("3 2\n", "Yes\n"), ("3 0\n", "No\n"), ("92 66\n", "Yes\n")])
