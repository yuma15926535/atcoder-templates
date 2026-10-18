# 最終照合の広い洗い出し（AtCoder_文字列問題一覧.txt の問題。作ったものを見ていないエージェントが解いたもの）
# agent_3 の担当問題（string_v2 最終照合）
# 実行: cd /home/claude/final/agent_3 && python3 /home/claude/final/run.py probs.py /mnt/user-data/uploads/claude/Main_base.cpp /home/claude/final/out /home/claude/ac-library [問題名…]
M2 = '/mnt/user-data/uploads/claude/math_v2/'

# ============================ C〜F ============================

# 中央が C の奇数長の部分文字列の数：各 'C' について min(左の数, 右の数) + 1
P('ABC458_C', [], r'''
void solv(){
    string S; cin >> S;
    ll n = S.size(), ans = 0;
    rep(i, n) if (S[i] == 'C') ans += min(i, n - 1 - i) + 1;
    cout << ans << '\n';
}
''', [("ABCCA\n", "5\n"), ("XYZ\n", "0\n"), ("SMBCPROGRAMMINGCONTEST\n", "11\n")])

# 同じ名前が前に何回出たか
P('ABC261_C', [], r'''
void solv(){
    int N; cin >> N;
    map<string, int> cnt;
    rep(i, N) {
        string s; cin >> s;
        int k = cnt[s]++;
        if (k == 0) cout << s << '\n';
        else cout << s << '(' << k << ")\n";
    }
}
''', [("5\nnewfile\nnewfile\nnewfolder\nnewfile\nnewfolder\n", "newfile\nnewfile(1)\nnewfolder\nnewfile(2)\nnewfolder(1)\n"),
      ("11\na\na\na\na\na\na\na\na\na\na\na\n", "a\na(1)\na(2)\na(3)\na(4)\na(5)\na(6)\na(7)\na(8)\na(9)\na(10)\n")])

# 部分列か：00_str_base の is_subsequence(t, s) と to_lower
P('ABC349_C', ['00_str_base'], r'''
void solv(){
    string S, T; cin >> S >> T;
    const string t = to_lower(T);
    if (is_subsequence(t, S) || (T[2] == 'X' && is_subsequence(t.substr(0, 2), S))) YES; else NO;
}
''', [("narita\nNRT\n", "Yes\n"), ("losangeles\nLAX\n", "Yes\n"), ("snuke\nRNG\n", "No\n")])

# 隣どうしが違う部分文字列の数 mod 998244353：隣が同じところで区切った区間ごとに len(len+1)/2
P('ABC456_C', [], r'''
void solv(){
    string S; cin >> S;
    ll n = S.size(), ans = 0, len = 0;
    rep(i, n) {
        len = (i > 0 && S[i] == S[i - 1]) ? 1 : len + 1;  // i で終わる「隣が違う」部分文字列の数
        ans += len;
    }
    cout << ans % 998244353 << '\n';
}
''', [("abbc\n", "6\n"),
      ("cabcabcbcaccacbcbcaabacbacaabccacbccbcacbacbacabcacabcaccaaaaabababcbabacaccabbcacbcbcbcababcbcbabca\n", "760\n")])

# A の挿入・削除：A を除いた骨格が同じなら、すき間ごとの A の数の差の和
P('ABC447_C', [], r'''
void solv(){
    string S, T; cin >> S >> T;
    auto gaps = [](const string& s) {
        string core; vector<ll> g(1, 0);
        for (char c : s) { if (c == 'A') g.back()++; else core += c, g.push_back(0); }
        return pair<string, vector<ll>>(core, g);
    };
    auto [cs, gs] = gaps(S);
    auto [ct, gt] = gaps(T);
    if (cs != ct) { cout << -1 << '\n'; return; }
    ll ans = 0;
    rep(i, sizei(gs)) ans += abs(gs[i] - gt[i]);
    cout << ans << '\n';
}
''', [("ABC\nBACA\n", "3\n"), ("ABC\nARC\n", "-1\n"), ("ABC\nABC\n", "0\n"), ("AAAWAZAAAABAAAU\nAWAAZABAAAAAUA\n", "9\n")])

# 文字の数を比べて、足りない分を @ で埋める（atcoder の文字だけ）：count_chars
P('ABC301_C', ['00_str_base'], r'''
void solv(){
    string S, T; cin >> S >> T;
    const Alphabet A = Alphabet::ascii();          // id(c) = 文字コード。'@' も入る
    const auto cs = count_chars(S, A), ct = count_chars(T, A);
    const string at = "atcoder";
    ll needS = 0, needT = 0;
    for (char c = 'a'; c <= 'z'; ++c) {
        const int d = cs[A.id(c)] - ct[A.id(c)];
        if (d == 0) continue;
        if (at.find(c) == string::npos) { NO; return; }
        if (d > 0) needT += d; else needS -= d;
    }
    if (needS <= cs[A.id('@')] && needT <= ct[A.id('@')]) YES; else NO;
}
''', [("ch@ku@ai\nchoku@@i\n", "Yes\n"), ("ch@kud@i\nakidu@ho\n", "Yes\n"), ("aoki\n@ok@\n", "No\n"), ("aa\nbb\n", "No\n")])

# 区間 [L, R] に T が現れるか：find_all ＋ lower_bound
P('ABC477_C', ['02_pattern_match'], r'''
void solv(){
    int Q; string S, T; cin >> Q >> S >> T;
    const vector<int> occ = find_all(S, T);       // T の開始位置（昇順）
    const int m = T.size();
    while (Q--) {
        int L, R; cin >> L >> R; --L;               // [L, R)
        auto it = lower_bound(occ.begin(), occ.end(), L);
        if (it != occ.end() && *it + m <= R) YES; else NO;
    }
}
''', [("3\nabcdabc\nbc\n2 6\n3 5\n3 7\n", "Yes\nNo\nYes\n"),
      ("10\nbbbaaabbbaaaaaabbbab\nbbb\n5 16\n5 6\n2 16\n6 7\n17 19\n3 7\n7 16\n7 20\n17 19\n2 6\n",
       "Yes\nNo\nYes\nNo\nNo\nNo\nYes\nYes\nNo\nNo\n")])

# 1 文字の書き換えと全体の大文字化・小文字化：最後の全体操作のときだけ to_lower / to_upper を全体にかける
P('ABC314_D', ['00_str_base'], r'''
void solv(){
    int N, Q; string S; cin >> N >> S >> Q;
    vector<tuple<int, int, char>> qs(Q);
    for (auto& [t, x, c] : qs) cin >> t >> x >> c;
    int g = -1;
    rep(i, Q) if (get<0>(qs[i]) != 1) g = i;
    rep(i, Q) {
        auto [t, x, c] = qs[i];
        if (t == 1) S[x - 1] = c;
        else if (i == g) S = (t == 2 ? to_lower(S) : to_upper(S));
    }
    cout << S << '\n';
}
''', [("7\nAtCoder\n5\n1 4 i\n3 0 a\n1 5 b\n2 0 a\n1 4 Y\n", "atcYber\n"),
      ("35\nTheQuickBrownFoxJumpsOverTheLazyDog\n10\n2 0 a\n1 19 G\n1 13 m\n1 2 E\n1 21 F\n2 0 a\n1 27 b\n3 0 a\n3 0 a\n1 15 i\n",
       "TEEQUICKBROWMFiXJUGPFOVERTBELAZYDOG\n")])

# 1 文字まで直せば回文になる部分文字列の数：中心ごとに RollingHash の lcp で「1 回目の不一致の先」まで飛ぶ（O(N log N)）
P('ABC468_D', ['01_rolling_hash'], r'''
void solv(){
    string S; cin >> S;
    const int n = S.size();
    const RollingHash H(S), HR(reversed(S));
    // 組 (b - k, a + k)（k = 0, 1, …）のうち、不一致が 1 個以下で続く組の数
    auto ext = [&](int a, int b) {
        const int maxk = min(b + 1, n - a);
        if (maxk <= 0) return 0;
        const int l1 = H.lcp(a, HR, n - 1 - b);   // s を b から左へ読んだもの = reversed の n-1-b から
        if (l1 + 1 >= maxk) return maxk;
        const int l2 = H.lcp(a + l1 + 1, HR, n - 1 - (b - l1 - 1));
        return min(maxk, l1 + 1 + l2);
    };
    ll ans = 0;
    rep(i, n) ans += 1 + ext(i + 1, i - 1);      // 奇数長（中心 i）
    rep(i, 1, n) ans += ext(i, i - 1);           // 偶数長（i-1 と i の間）
    cout << ans << '\n';
}
''', [("ababa\n", "13\n"), ("atcoder\n", "18\n"), ("abccbacbacb\n", "40\n")])

# 部分文字列を 1 つ左に回して辞書順最小：最初の下がる位置の文字を、それより大きい文字の手前まで運ぶ
P('ABC409_D', [], r'''
void solv(){
    int T; cin >> T;
    while (T--) {
        int N; string S; cin >> N >> S;
        int i = 0;
        while (i + 1 < N && S[i] <= S[i + 1]) ++i;
        if (i + 1 >= N) { cout << S << '\n'; continue; }
        const char c = S[i];
        int j = i + 1;
        while (j < N && S[j] <= c) ++j;
        cout << S.substr(0, i) + S.substr(i + 1, j - i - 1) + c + S.substr(j) << '\n';
    }
}
''', [("3\n7\natcoder\n1\nx\n5\nsnuke\n", "acodert\nx\nnsuke\n")])

# 隣が等しい所がちょうど 1 か所：0101… の型の接頭辞の費用 ＋ 反対の型の接尾辞の費用
P('ABC346_D', [], r'''
void solv(){
    int N; string S; cin >> N >> S;
    vector<ll> C(N); cin >> C;
    vector<array<ll, 2>> pre(N + 1, {0, 0}), suf(N + 1, {0, 0});
    rep(i, N) rep(p, 2) pre[i + 1][p] = pre[i][p] + ((S[i] - '0') != ((i & 1) ^ p) ? C[i] : 0);
    rrep(i, N) rep(p, 2) suf[i][p] = suf[i + 1][p] + ((S[i] - '0') != ((i & 1) ^ p) ? C[i] : 0);
    ll ans = LINF;
    rep(i, 1, N) rep(p, 2) chmin(ans, pre[i][p] + suf[i][p ^ 1]);   // i-1 と i が等しい
    cout << ans << '\n';
}
''', [("5\n00011\n3 9 2 6 4\n", "7\n"), ("4\n1001\n1 2 3 4\n", "0\n"),
      ("11\n11111100111\n512298012 821282085 543342199 868532399 690830957 973970164 928915367 954764623 923012648 540375785 925723427\n", "2286846953\n")])

# T を部分列に含まない部分文字列の数：各開始位置 l で NextTable::match_end（含む最短の終わり）
P('ABC452_D', ['06_subsequence'], r'''
void solv(){
    string S, T; cin >> S >> T;
    const int n = S.size();
    const NextTable nx(S);
    ll ans = 0;
    rep(l, n) {
        const int e = nx.match_end(T, l);          // S[l, e) が T を含む最短（なければ -1）
        ans += (e < 0 ? n - l : e - 1 - l);
    }
    cout << ans << '\n';
}
''', [("abrakadabra\naba\n", "51\n"), ("aaaaa\na\n", "0\n"), ("rdddrdtdcdrrdcredctdordoeecrotet\ndcre\n", "263\n")])

# a が A 個・b が B 個の並べ替えの辞書順 K 番目：math_v2 04_permutation の multiset_perm_unrank
P('ABC202_D', [M2 + '04_permutation'], r'''
void solv(){
    int A, B; ll K; cin >> A >> B >> K;
    cout << multiset_perm_unrank(string(A, 'a') + string(B, 'b'), (Int128)(K - 1)) << '\n';
}
''', [("2 2 4\n", "baab\n"), ("30 30 118264581564861424\n", "bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa\n")])

# 全部の部分文字列の値の和（巨大）：10^p の係数を累積和で作り、math_v2 01_digits の carry_to_string で繰り上げる
P('ABC379_E', [M2 + '01_digits'], r'''
void solv(){
    int N; string S; cin >> N >> S;
    vector<ll> pre(N + 1, 0);                         // pre[i] = Σ_{j<i} S[j] × (j + 1)
    rep(i, N) pre[i + 1] = pre[i] + (ll)(S[i] - '0') * (i + 1);
    vector<ll> c(N);                                  // c[p] = 10^p の係数
    rep(p, N) c[p] = pre[N - p];
    cout << carry_to_string(c) << '\n';
}
''', [("3\n379\n", "514\n"), ("30\n314159265358979323846264338327\n", "369673254065355789035427227741\n")])

# S_n = S_{n-1} + S_{n-2} の区間の文字の数：長さが 10^18 を超えるまで作り、接頭辞の数を再帰（ループ）で。土台は CharPrefix
P('ABC450_E', ['00_str_base'], r'''
void solv(){
    string X, Y; cin >> X >> Y;
    const CharPrefix cx(X), cy(Y);
    vector<ll> len = {0, (ll)X.size(), (ll)Y.size()};
    vector<vector<ll>> tot = {vector<ll>(26), vector<ll>(26), vector<ll>(26)};
    rep(c, 26) tot[1][c] = cx.total((char)('a' + c)), tot[2][c] = cy.total((char)('a' + c));
    const ll LIM = 1000000000000000000LL;
    while (len.back() < LIM) {
        const int k = len.size();
        len.push_back(len[k - 1] + len[k - 2]);
        vector<ll> t(26);
        rep(c, 26) t[c] = tot[k - 1][c] + tot[k - 2][c];
        tot.push_back(t);
    }
    const int top = max<int>(2, (int)len.size() - 1);   // S_top は S_{10^18} の接頭辞（長さ ≥ 10^18）
    auto pre = [&](ll r, char ch) {                       // S_top の先頭 r 文字の ch の数
        int k = top; ll res = 0; const int c = ch - 'a';
        while (k > 2) {
            if (r <= len[k - 1]) --k;
            else res += tot[k - 1][c], r -= len[k - 1], k -= 2;
        }
        return res + (k == 1 ? cx.count(0, (int)r, ch) : cy.count(0, (int)r, ch));
    };
    int Q; cin >> Q;
    while (Q--) {
        ll L, R; char C; cin >> L >> R >> C;
        cout << pre(R, C) - pre(L - 1, C) << '\n';
    }
}
''', [("a\nb\n6\n2 7 a\n1 3 b\n3 7 b\n1 9 c\n1 1000000000000000000 b\n1000000000000000000 1000000000000000000 a\n",
       "3\n2\n3\n0\n618033988749894848\n1\n")])

# 11/22 の部分列の最長（区間クエリ）：CharPrefix で区間の 1・2・/ の数、'/' の位置で二分探索
P('ABC381_E', ['00_str_base'], r'''
void solv(){
    int N, Q; string S; cin >> N >> Q >> S;
    const CharPrefix cp(S, Alphabet::of("12/"));
    while (Q--) {
        int L, R; cin >> L >> R; --L;                       // [L, R)
        int lo = cp.count(0, L, '/'), hi = cp.count(0, R, '/');   // 区間の '/' は lo〜hi-1 番目
        if (lo == hi) { cout << 0 << '\n'; continue; }
        auto ones = [&](int k) { return cp.count(L, cp.kth('/', k), '1'); };
        auto twos = [&](int k) { return cp.count(cp.kth('/', k) + 1, R, '2'); };
        int a = lo, b = hi - 1;                             // ones(k) ≥ twos(k) となる最初の k
        while (a < b) { int m = (a + b) / 2; if (ones(m) >= twos(m)) b = m; else a = m + 1; }
        ll best = min(ones(a), twos(a));
        if (a > lo) chmax(best, min(ones(a - 1), twos(a - 1)));
        cout << 2 * best + 1 << '\n';
    }
}
''', [("12 5\n111/212/1122\n1 7\n9 12\n3 6\n4 10\n1 12\n", "5\n0\n3\n1\n7\n")])

# ============================ B ============================

# 4 つが全部違うか
P('ABC211_B', [], r'''
void solv(){
    set<string> st;
    rep(i, 4) { string s; cin >> s; st.insert(s); }
    if (st.size() == 4) YES; else NO;
}
''', [("3B\nHR\n2B\nH\n", "Yes\n"), ("2B\n3B\nHR\n3B\n", "No\n")])

# 区間の反転
P('ABC233_B', [], r'''
void solv(){
    int L, R; string S; cin >> L >> R >> S;
    reverse(S.begin() + (L - 1), S.begin() + R);
    cout << S << '\n';
}
''', [("3 7\nabcdefgh\n", "abgfedch\n"), ("1 7\nreviver\n", "reviver\n"), ("4 13\nmerrychristmas\n", "meramtsirhcyrs\n")])

# S が T の接頭辞か：common_prefix_length
P('ABC268_B', ['00_str_base'], r'''
void solv(){
    string S, T; cin >> S >> T;
    if (common_prefix_length(S, T) == sizei(S)) YES; else NO;
}
''', [("atco\natcoder\n", "Yes\n"), ("code\natcoder\n", "No\n"), ("abc\nabc\n", "Yes\n"), ("aaaa\naa\n", "No\n")])

# 最初に違う位置：common_prefix_length + 1
P('ABC374_B', ['00_str_base'], r'''
void solv(){
    string S, T; cin >> S >> T;
    cout << (S == T ? 0 : common_prefix_length(S, T) + 1) << '\n';
}
''', [("abcde\nabedc\n", "3\n"), ("abcde\nabcdefg\n", "6\n"), ("keyence\nkeyence\n", "0\n")])

# 中央ぞろえ
P('ABC445_B', [], r'''
void solv(){
    int N; cin >> N;
    vector<string> S(N); cin >> S;
    size_t m = 0;
    for (auto& s : S) m = max(m, s.size());
    for (auto& s : S) { string p((m - s.size()) / 2, '.'); cout << p + s + p << '\n'; }
}
''', [("4\napple\nblueberry\ncoconut\ndragonfruit\n", "...apple...\n.blueberry.\n..coconut..\ndragonfruit\n"),
      ("6\nabc\nd\nefghi\njkl\nmnopq\nr\n", ".abc.\n..d..\nefghi\n.jkl.\nmnopq\n..r..\n")])

# oxx の繰り返しの部分文字列か：find_all_cyclic
P('ABC230_B', ['02_pattern_match'], r'''
void solv(){
    string S; cin >> S;
    if (!find_all_cyclic("oxx", S).empty()) YES; else NO;
}
''', [("xoxxoxxo\n", "Yes\n"), ("xxoxxoxo\n", "No\n"), ("ox\n", "Yes\n")])

# 連続する O の長さ / K の和：run_length
P('ABC379_B', ['00_str_base'], r'''
void solv(){
    int N, K; string S; cin >> N >> K >> S;
    ll ans = 0;
    for (const auto& r : run_length(S)) if (r.value == 'O') ans += r.len / K;
    cout << ans << '\n';
}
''', [("7 3\nOOXOOOO\n", "1\n"), ("12 2\nOXXOOOXOOOOX\n", "3\n"), ("22 5\nXXOOOOOOOOXXOOOOOXXXXX\n", "2\n")])

# 上位 K 人を辞書順に
P('ABC288_B', [], r'''
void solv(){
    int N, K; cin >> N >> K;
    vector<string> S(N); cin >> S;
    sort(S.begin(), S.begin() + K);
    rep(i, K) cout << S[i] << '\n';
}
''', [("5 3\nabc\naaaaa\nxyz\na\ndef\n", "aaaaa\nabc\nxyz\n"), ("4 4\nz\nzyx\nzzz\nrbg\n", "rbg\nz\nzyx\nzzz\n"), ("3 1\nabc\narc\nagc\n", "abc\n")])

# 大文字・小文字がある・全部違う：to_lower / to_upper と count_chars(S, Alphabet::letters())
P('ABC249_B', ['00_str_base'], r'''
void solv(){
    string S; cin >> S;
    const auto cnt = count_chars(S, Alphabet::letters());
    const bool distinct = *max_element(cnt.begin(), cnt.end()) <= 1;
    if (to_lower(S) != S && to_upper(S) != S && distinct) YES; else NO;
}
''', [("AtCoder\n", "Yes\n"), ("Aa\n", "Yes\n"), ("atcoder\n", "No\n"), ("Perfect\n", "No\n")])

# 大文字・小文字を区別しない最頻：to_lower
P('ABC471_B', ['00_str_base'], r'''
void solv(){
    int N; cin >> N;
    map<string, int> cnt;
    int ans = 0;
    rep(i, N) { string s; cin >> s; chmax(ans, ++cnt[to_lower(s)]); }
    cout << ans << '\n';
}
''', [("7\nARC\nabc\nahc\nABC\nbeginner\nAbC\nahc\n", "3\n"), ("10\nx\nx\nx\nx\nx\nx\nx\nx\nx\nx\n", "10\n")])

# 1122 文字列：run_length がすべて長さ 2 で、文字が全部違う
P('ABC381_B', ['00_str_base'], r'''
void solv(){
    string S; cin >> S;
    const auto rl = run_length(S);
    set<char> st;
    bool ok = true;
    for (const auto& r : rl) ok &= r.len == 2, st.insert(r.value);
    if (ok && st.size() == rl.size()) YES; else NO;
}
''', [("aabbcc\n", "Yes\n"), ("aab\n", "No\n"), ("zzzzzz\n", "No\n")])

# キーボードの移動距離：Alphabet::of(S).id(c) が c の位置
P('ABC373_B', ['00_str_base'], r'''
void solv(){
    string S; cin >> S;
    const Alphabet K = Alphabet::of(S);
    ll ans = 0;
    for (char c = 'A'; c < 'Z'; ++c) ans += abs(K.id(c + 1) - K.id(c));
    cout << ans << '\n';
}
''', [("ABCDEFGHIJKLMNOPQRSTUVWXYZ\n", "25\n"), ("MGJYIZDKSBHPVENFLQURTCWOAX\n", "223\n")])

# # で区切った区間ごとに 1 個の o：split ＋ join（答えは複数あるので check で確かめる）
def chk_416(inp, out, exp):
    S = inp.split()[0]
    o = out.split()
    if len(o) != 1: return False
    T, E = o[0], exp.split()[0]
    if len(T) != len(S): return False
    for a, b in zip(S, T):
        if (a == '#') != (b == '#') or b not in '.#o': return False
    return all(seg.count('o') <= 1 for seg in T.split('#')) and T.count('o') == E.count('o')
P('ABC416_B', ['00_str_base'], r'''
void solv(){
    string S; cin >> S;
    auto parts = split(S, '#');                 // 空の部分も残る（# の位置を保つ）
    for (auto& p : parts) if (!p.empty()) p[0] = 'o';
    cout << join(parts, "#") << '\n';
}
''', [("#..#.\n", "#o.#o\n"), ("#\n", "#\n"), (".....\n", "..o..\n"), ("...#..#.##.#.\n", "o..#.o#o##o#o\n")], check=chk_416)

# i と o を交互にするための最小の挿入
P('ABC397_B', [], r'''
void solv(){
    string S; cin >> S;
    char e = 'i'; int ans = 0;
    for (char c : S) {
        if (c != e) ++ans, e = (e == 'i' ? 'o' : 'i');   // e を挿入
        e = (e == 'i' ? 'o' : 'i');                       // c を使う
    }
    if (e == 'o') ++ans;
    cout << ans << '\n';
}
''', [("ioi\n", "1\n"), ("iioo\n", "2\n"), ("io\n", "0\n")])

# S_i + S_j が回文になる組：is_palindrome
P('ABC307_B', ['00_str_base'], r'''
void solv(){
    int N; cin >> N;
    vector<string> S(N); cin >> S;
    rep(i, N) rep(j, N) if (i != j && is_palindrome(S[i] + S[j])) { YES; return; }
    NO;
}
''', [("5\nab\nccef\nda\na\nfe\n", "Yes\n"), ("3\na\nb\naba\n", "No\n"),
      ("2\naaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa\naaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa\n", "Yes\n")])

# 幅 w で折り返して c 列目を縦に読む：全探索
P('ABC360_B', [], r'''
void solv(){
    string S, T; cin >> S >> T;
    const int n = S.size();
    rep(w, 1, n) rep(c, w) {
        string t;
        for (int k = c; k < n; k += w) t += S[k];
        if (t == T) { YES; return; }
    }
    NO;
}
''', [("atcoder toe\n", "Yes\n"), ("beginner r\n", "No\n"), ("verticalreading agh\n", "No\n")])

# 数字を +1 (9→0) して部分文字列にする最小回数：ずらしごとに和
P('ABC438_B', [], r'''
void solv(){
    int N, M; string S, T; cin >> N >> M >> S >> T;
    int ans = INF;
    rep(i, N - M + 1) {
        int s = 0;
        rep(j, M) s += ((S[i + j] - T[j]) % 10 + 10) % 10;
        chmin(ans, s);
    }
    cout << ans << '\n';
}
''', [("4 2\n2025\n91\n", "2\n"), ("3 2\n438\n38\n", "0\n"), ("5 5\n00000\n11111\n", "45\n"), ("8 3\n20251227\n438\n", "13\n")])

# 長さ K の部分文字列の最頻とその一覧（辞書順）：map で数えて join で出力
P('ABC428_B', ['00_str_base'], r'''
void solv(){
    int N, K; string S; cin >> N >> K >> S;
    map<string, int> cnt;
    rep(i, N - K + 1) cnt[S.substr(i, K)]++;
    int x = 0;
    for (auto& [s, c] : cnt) chmax(x, c);
    vector<string> ans;
    for (auto& [s, c] : cnt) if (c == x) ans.push_back(s);
    cout << x << '\n' << join(ans) << '\n';
}
''', [("9 3\novowowovo\n", "2\novo owo\n"), ("9 1\novowowovo\n", "5\no\n"), ("35 3\nthequickbrownfoxjumpsoverthelazydog\n", "2\nthe\n")])

# あだ名：自分の姓か名が、ほかの人の姓・名のどれとも一致しないか
P('ABC247_B', [], r'''
void solv(){
    int N; cin >> N;
    vector<string> s(N), t(N);
    rep(i, N) cin >> s[i] >> t[i];
    rep(i, N) {
        bool okS = true, okT = true;
        rep(j, N) if (j != i) {
            if (s[i] == s[j] || s[i] == t[j]) okS = false;
            if (t[i] == s[j] || t[i] == t[j]) okT = false;
        }
        if (!okS && !okT) { NO; return; }
    }
    YES;
}
''', [("3\ntanaka taro\ntanaka jiro\nsuzuki hanako\n", "Yes\n"), ("3\naaa bbb\nxxx aaa\nbbb yyy\n", "No\n"),
      ("2\ntanaka taro\ntanaka taro\n", "No\n"), ("3\ntakahashi chokudai\naoki kensho\nsnu ke\n", "Yes\n")])
