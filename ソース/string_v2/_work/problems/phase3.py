# Phase 3（★ 11 本）の代表的な問題。入出力例は AtCoder の問題ページで確かめたもの
M = '/mnt/user-data/uploads/claude/math_v2/00_modint'  # ModInt107 / ModInt998 を使うとき

# ---------------- 答えが複数あるときの確かめ方 ----------------
def check_abc284f(inp, out, exp):
    n, t = inp.split()[0], inp.split()[1]
    n = int(n)
    o = out.split()
    if exp.split() == ['-1']: return o == ['-1']
    if len(o) != 2: return False
    s, i = o[0], int(o[1])
    return len(s) == n and 0 <= i <= n and s[:i] + s[::-1] + s[i:] == t

def check_lcs(inp, out, exp):
    a, b = (inp.split('\n') + ['', ''])[:2]
    o = out.strip()
    def sub(x, y):
        it = iter(y)
        return all(c in it for c in x)
    return len(o) == len(exp.strip()) and sub(o, a) and sub(o, b)

def check_abc280h(inp, out, exp):
    tok = inp.split()
    n = int(tok[0]); ss = tok[1:1 + n]
    def get(k, l, r): return ss[int(k) - 1][int(l) - 1:int(r)]
    o = out.split(); e = exp.split()
    if len(o) != len(e): return False
    for p in range(0, len(e), 3):
        if get(*o[p:p + 3]) != get(*e[p:p + 3]): return False
    return True

# ================= 01_rolling_hash =================
P('ABC141_E_RH', ['01_rolling_hash'], r'''
void solv(){
    int N; string S; cin >> N >> S;
    RollingHash H(S);
    // 長さ len で重ならずに 2 回現れるか（初めて見た位置を覚える）
    auto ok = [&](int len) {
        unordered_map<StrHash, int> first;
        for (int i = 0; i + len <= N; ++i) {
            auto [it, inserted] = first.emplace(H.get(i, i + len), i);
            if (!inserted && i - it->second >= len) return true;
        }
        return false;
    };
    int lo = 0, hi = N / 2;
    while (lo < hi) { int mid = (lo + hi + 1) / 2; if (ok(mid)) lo = mid; else hi = mid - 1; }
    cout << lo << '\n';
}''', [("5\nababa\n", "2\n"), ("2\nxy\n", "0\n"), ("13\nstrangeorange\n", "5\n")])

P('ABC284_F', ['01_rolling_hash'], r'''
void solv(){
    int N; string T; cin >> N >> T;
    const string R = reversed(T);
    RollingHash H(T), HR(R);
    // T = S[0, i) + rev(S) + S[i, N) で S = T[0, i) + T[i+N, 2N)
    for (int i = 0; i <= N; ++i) {
        if (H.get(i, N) == HR.get(0, N - i) && H.get(N, N + i) == HR.get(2 * N - i, 2 * N)) {
            cout << T.substr(0, i) + T.substr(i + N) << '\n' << i << '\n';
            return;
        }
    }
    cout << -1 << '\n';
}''', [("3\nabcbac\n", "abc\n2\n"), ("4\nabababab\n", "abab\n1\n"), ("3\nagccga\n", "cga\n0\n"), ("4\natcodeer\n", "-1\n")], check_abc284f)

P('TYP047', ['01_rolling_hash'], r'''
void solv(){
    int N; string S, T; cin >> N >> S >> T;
    auto id = [](char c) { return c == 'R' ? 0 : c == 'G' ? 1 : 2; };
    vector<int> a(N);
    rep(i, N) a[i] = id(S[i]);
    RollingHash HS(a);
    // マス (i, j) の色が c ⇔ s_i + t_j + c ≡ 0 (mod 3)
    vector<RollingHash> HT;
    rep(c, 3) {
        vector<int> b(N);
        rep(j, N) b[j] = ((-id(T[j]) - (int)c) % 3 + 3) % 3;
        HT.emplace_back(b);
    }
    ll ans = 0;
    for (int k = -(N - 1); k <= N - 1; ++k) {
        bool any = false;
        rep(c, 3) {
            const auto& B = HT[c];
            if (k >= 0) any |= HS.get(0, N - k) == B.get(k, N);
            else any |= HS.get(-k, N) == B.get(0, N + k);
        }
        ans += any;
    }
    cout << ans << '\n';
}''', [("5\nRGBGB\nGRGRB\n", "6\n"), ("3\nRRR\nBBB\n", "5\n"), ("10\nBGGGRBBGRG\nRGBBRGRGGG\n", "4\n")])

# ================= 02_pattern_match =================
P('ABC257_G', ['02_pattern_match'], r'''
void solv(){
    string S, T; cin >> S >> T;
    const auto L = match_lengths(T, S);  // T[i..) と S の最長共通接頭辞
    const int n = sizei(T);
    // 0 から n まで、i から i+1〜i+L[i] に飛べるときの最小の回数（貪欲）
    int cnt = 0, cur_end = 0, far = 0;
    for (int i = 0; i < n; ++i) {
        if (i > cur_end) { cout << -1 << '\n'; return; }
        chmax(far, i + L[i]);
        if (i == cur_end) {
            if (far == cur_end) { cout << -1 << '\n'; return; }
            ++cnt, cur_end = far;
        }
    }
    cout << (cur_end >= n ? cnt : -1) << '\n';
}''', [("aba\nababaab\n", "3\n"), ("atcoder\nac\n", "-1\n")])

P('ABC150_F', ['02_pattern_match'], r'''
void solv(){
    int N; cin >> N;
    vector<ll> a(N), b(N); cin >> a >> b;
    vector<ll> da(N), db(N);
    rep(i, N) da[i] = a[i] ^ a[(i + 1) % N], db[i] = b[i] ^ b[(i + 1) % N];
    for (int k : rotation_offsets(da, db)) cout << k << ' ' << (a[k] ^ b[0]) << '\n';
}''', [("3\n0 2 1\n1 2 3\n", "1 3\n"), ("5\n0 0 0 0 0\n2 2 2 2 2\n", "0 2\n1 2\n2 2\n3 2\n4 2\n"),
       ("6\n0 1 3 7 6 4\n1 5 4 6 2 3\n", "2 2\n5 5\n"), ("2\n1 2\n0 0\n", "")])

P('ABC135_F', ['02_pattern_match'], r'''
void solv(){
    string s, t; cin >> s >> t;
    const int n = sizei(s), m = sizei(t);
    vector<char> occ(n, 0);
    for (int i : find_all_cyclic(s, t)) occ[i] = 1;
    // i → (i + m) mod n（occ[i] のとき）。続けて現れる最大の回数。輪があれば無限
    vector<int> len(n, -1), state(n, 0);  // state: 0 未訪問, 1 たどり中, 2 済み
    int ans = 0;
    for (int s0 = 0; s0 < n; ++s0) {
        if (state[s0]) continue;
        vector<int> path;
        int v = s0;
        while (state[v] == 0 && occ[v]) state[v] = 1, path.push_back(v), v = (v + m) % n;
        if (occ[v] && state[v] == 1) { cout << -1 << '\n'; return; }
        int base = occ[v] ? len[v] : 0;
        if (!occ[v]) state[v] = 2, len[v] = 0;
        for (int k = sizei(path) - 1; k >= 0; --k) len[path[k]] = ++base, state[path[k]] = 2, chmax(ans, base);
    }
    cout << ans << '\n';
}''', [("abcabab\nab\n", "3\n"), ("aa\naaaaaaa\n", "-1\n"), ("aba\nbaaab\n", "0\n")])

P('ABC430_E', ['02_pattern_match'], r'''
void solv(){
    int T; cin >> T;
    while (T--) {
        string A, B; cin >> A >> B;
        const auto ks = rotation_offsets(A, B);
        cout << (ks.empty() ? -1 : ks[0]) << '\n';
    }
}''', [("5\n1010001\n1000110\n000\n111\n01010\n01010\n0101\n0011\n100001101110000001010110110001\n101100011000011011100000010101\n", "2\n-1\n0\n-1\n22\n")])

P('ABC343_G', ['02_pattern_match'], r'''
void solv(){
    int N; cin >> N;
    vector<string> v(N); cin >> v;
    sort(nall(v)); v.erase(unique(nall(v)), v.end());
    vector<string> w;  // ほかに含まれるものを除く
    rep(i, sizei(v)) {
        bool inside = false;
        rep(j, sizei(v)) if (i != j && count_occurrences(v[j], v[i]) > 0) inside = true;
        if (!inside) w.push_back(v[i]);
    }
    const int n = sizei(w);
    vector<vector<int>> ov(n, vector<int>(n));
    rep(i, n) rep(j, n) if (i != j) ov[i][j] = max_overlap(w[i], w[j]);
    vector<vector<ll>> dp(1 << n, vector<ll>(n, LINF));
    rep(i, n) dp[1 << i][i] = sizei(w[i]);
    rep(S, 1 << n) rep(i, n) if (dp[S][i] < LINF) rep(j, n) if (!(S >> j & 1)) chmin(dp[S | 1 << j][j], dp[S][i] + sizei(w[j]) - ov[i][j]);
    cout << *min_element(nall(dp[(1 << n) - 1])) << '\n';
}''', [("3\nsnuke\nkensho\nuk\n", "9\n"), ("3\nabc\nabc\narc\n", "6\n"), ("6\ncmcmrcc\nrmrrrmr\nmrccm\nmmcr\nrmmrmrcc\nccmcrcmcm\n", "27\n")])

# ================= 03_suffix_array =================
P('ABC141_E_SA', ['03_suffix_array'], r'''
void solv(){
    int N; string S; cin >> N >> S;
    cout << SuffixArray(S).longest_repeated_non_overlapping().len << '\n';
}''', [("5\nababa\n", "2\n"), ("2\nxy\n", "0\n"), ("13\nstrangeorange\n", "5\n")])

P('ABC213_F', ['03_suffix_array'], r'''
void solv(){
    int N; string S; cin >> N >> S;
    for (ll x : SuffixArray(S).sum_lcp_with_all()) cout << x << '\n';
}''', [("3\nabb\n", "3\n3\n2\n"), ("11\nmississippi\n", "11\n16\n14\n12\n13\n11\n9\n7\n4\n3\n4\n")])

P('ABC272_F', ['03_suffix_array'], r'''
void solv(){
    int N; string S, T; cin >> N >> S >> T;
    // S S '$' T T '~' の接尾辞配列。同じ回転なら S の方が前に来る（'$' はどの文字より小さく、'~' は大きい）
    const string X = S + S + '$' + T + T + '~';
    ll ans = 0, seen = 0;
    for (int p : make_suffix_array(X)) {
        if (p < N) ++seen;
        else if (2 * N + 1 <= p && p < 3 * N + 1) ans += seen;
    }
    cout << ans << '\n';
}''', [("3\nadb\ncab\n", "4\n"), ("10\nwsiuhwijsl\npwqoketvun\n", "56\n")])

P('ABC280_Ex', ['03_suffix_array'], r'''
void solv(){
    int N; cin >> N;
    vector<string> S(N); cin >> S;
    const auto SA = SuffixArray::of_strings(S);
    int Q; cin >> Q;
    while (Q--) {
        ll x; cin >> x;
        const StrSpan sp = SA.kth_substring_with_dup(x - 1);
        cout << SA.owner(sp.l) + 1 << ' ' << SA.offset(sp.l) + 1 << ' ' << SA.offset(sp.l) + sp.len() << '\n';
    }
}''', [("2\nabab\ncab\n2\n5 14\n", "1 3 4\n2 1 1\n"), ("3\na\na\nba\n2\n1 2\n", "1 1 1\n1 1 1\n"),
       ("10\ngxgpuamkx\nszhkbpphykin\nezplvfja\nmopodotkrj\nrimlvumuar\nnexcfyce\neurgvjyos\ndhvuyfvt\nnrdyluacvra\nggwnpnzij\n6\n74 268 310 380 455 489\n",
        "3 1 2\n4 4 5\n4 3 7\n9 6 6\n6 6 6\n2 2 12\n")], check_abc280h)

# ================= 04_palindrome =================
P('ABC237_C', ['04_palindrome'], r'''
void solv(){
    string S; cin >> S;
    const int n = sizei(S);
    int p = 0, q = 0;
    while (p < n && S[p] == 'a') ++p;
    if (p == n) { YES; return; }
    while (S[n - 1 - q] == 'a') ++q;
    if (q >= p && Palindromes(S).is_palindrome(p, n - q)) YES; else NO;
}''', [("kasaka\n", "Yes\n"), ("atcoder\n", "No\n"), ("php\n", "Yes\n")])

P('ABC320_B', ['04_palindrome'], r'''
void solv(){
    string S; cin >> S;
    cout << Palindromes(S).longest().len() << '\n';
}''', [("TOYOTA\n", "5\n"), ("ABCDEFG\n", "1\n"), ("AAAAAAAAAA\n", "10\n")])

P('ABC398_F', ['04_palindrome'], r'''
void solv(){
    string S; cin >> S;
    cout << min_palindrome_with_prefix(S) << '\n';
}''', [("ABC\n", "ABCBA\n"), ("Z\n", "Z\n"), ("TREE\n", "TREERT\n")])

# ================= 05_lex_order =================
P('ABC225_F', ['05_lex_order'], r'''
void solv(){
    int N, K; cin >> N >> K;
    vector<string> S(N); cin >> S;
    vector<string> v;
    for (int i : concat_order(S)) v.push_back(S[i]);
    // 並べた順で、後ろから「i 以降から k 枚」の辞書順最小
    vector<vector<string>> dp(N + 1, vector<string>(K + 1, "~"));
    rep(i, N + 1) dp[i][0] = "";
    for (int i = N - 1; i >= 0; --i)
        for (int k = 1; k <= K; ++k) {
            dp[i][k] = dp[i + 1][k];
            if (dp[i + 1][k - 1] != "~") chmin(dp[i][k], v[i] + dp[i + 1][k - 1]);
        }
    cout << dp[0][K] << '\n';
}''', [("4 3\node\nzaaa\nr\natc\n", "atcoder\n"), ("5 2\nz\nz\nzzz\nz\nzzzzzz\n", "zz\n")])

P('ABC219_C', ['05_lex_order'], r'''
void solv(){
    string X; int N; cin >> X >> N;
    vector<string> S(N); cin >> S;
    sort(nall(S), lex_less_by(X));
    for (auto& s : S) cout << s << '\n';
}''', [("bacdefghijklmnopqrstuvwxzy\n4\nabx\nbzz\nbzy\ncaa\n", "bzz\nbzy\nabx\ncaa\n"), ("zyxwvutsrqponmlkjihgfedcba\n5\na\nab\nabc\nac\nb\n", "b\na\nac\nab\nabc\n")])

P('ABC171_C', ['05_lex_order'], r'''
void solv(){
    ll N; cin >> N;
    cout << bijective_to_string(N) << '\n';
}''', [("2\n", "b\n"), ("27\n", "aa\n"), ("123456789\n", "jjddja\n")])

P('ABC285_C', ['05_lex_order'], r'''
void solv(){
    string S; cin >> S;
    cout << bijective_from_string(S, Alphabet::upper()) << '\n';
}''', [("AB\n", "28\n"), ("C\n", "3\n"), ("BRUTMHYHIIZP\n", "10000000000000000\n")])

P('ABC223_B', ['05_lex_order'], r'''
void solv(){
    string S; cin >> S;
    const int a = min_rotation(S), b = max_rotation(S);
    cout << S.substr(a) + S.substr(0, a) << '\n' << S.substr(b) + S.substr(0, b) << '\n';
}''', [("aaba\n", "aaab\nbaaa\n"), ("z\n", "z\nz\n"), ("abracadabra\n", "aabracadabr\nracadabraab\n")])

# ================= 06_subsequence =================
P('TYP006', ['06_subsequence'], r'''
void solv(){
    int N, K; string S; cin >> N >> K >> S;
    cout << lex_smallest_subsequence(S, K) << '\n';
}''', [("7 3\natcoder\n", "acd\n"), ("14 5\nkittyonyourlap\n", "inlap\n")])

P('TYP008', ['06_subsequence', M], r'''
void solv(){
    int N; string S; cin >> N >> S;
    cout << count_subsequence_occurrences<ModInt107>(S, "atcoder") << '\n';
}''', [("10\nattcordeer\n", "4\n"), ("41\nbtwogablwetwoiehocghiewobadegwhoihegnldir\n", "2\n"),
       ("140\n" + ''.join(c * 20 for c in "atcoder") + "\n", "279999993\n")])

P('ABC211_C', ['06_subsequence', M], r'''
void solv(){
    string S; cin >> S;
    cout << count_subsequence_occurrences<ModInt107>(S, "chokudai") << '\n';
}''', [("chchokudai\n", "3\n"), ("atcoderrr\n", "0\n"), ("chokudaichokudaichokudai\n", "45\n")])

P('ABC214_F', ['06_subsequence', M], r'''
void solv(){
    string S; cin >> S;
    cout << count_distinct_subsequences<ModInt107>(S, false, 2) << '\n';
}''', [("abc\n", "4\n"), ("aa\n", "1\n"), ("acba\n", "6\n"), ("chokudai\n", "54\n")])

P('ABC230_F', ['06_subsequence', M], r'''
void solv(){
    int N; cin >> N;
    vector<ll> A(N); cin >> A;
    // 累積和 P_1〜P_{N-1} の異なる部分列（空も）と 1 対 1
    vector<ll> P;
    ll s = 0;
    rep(i, N - 1) s += A[i], P.push_back(s);
    cout << count_distinct_subsequences<ModInt998>(P, true) << '\n';
}''', [("3\n1 -1 1\n", "4\n"), ("10\n377914575 -275478149 0 -444175904 719654053 -254224494 -123690081 377914575 -254224494 -21253655\n", "321\n")])

P('ARC081_E', ['06_subsequence'], r'''
void solv(){
    string A; cin >> A;
    cout << shortest_non_subsequence(A) << '\n';
}''', [("atcoderregularcontest\n", "b\n"), ("abcdefghijklmnopqrstuvwxyz\n", "aa\n")])

P('ABC324_E', ['06_subsequence'], r'''
void solv(){
    int N; string T; cin >> N >> T;
    const int m = sizei(T);
    vector<int> pre(N), suf(N);
    rep(i, N) {
        string S; cin >> S;
        pre[i] = prefix_match_lengths(S, T).back();  // T の先頭何文字までが S の部分列か
        suf[i] = suffix_match_lengths(S, T)[0];      // T の末尾何文字までが S の部分列か
    }
    sort(nall(suf));
    ll ans = 0;
    rep(i, N) ans += suf.end() - lower_bound(nall(suf), m - pre[i]);
    cout << ans << '\n';
}''', [("3 bac\nabba\nbcb\naaca\n", "3\n"), ("5 xx\nx\nx\nx\nx\nx\n", "25\n"), ("1 y\nx\n", "0\n"),
       ("10 ms\nmkgn\nm\nhlms\nvmsle\nmxsm\nnnzdhi\numsavxlb\nffnsybomr\nyvmm\nnaouel\n", "68\n")])

# ================= 06_lcs_edit =================
P('EDPC_F', ['06_lcs_edit'], r'''
void solv(){
    string s, t; cin >> s >> t;
    cout << lcs_restore(s, t) << '\n';
}''', [("axyb\nabyxb\n", "axb\n"), ("aa\nxayaz\n", "aa\n"), ("a\nz\n", "\n"), ("abracadabra\navadakedavra\n", "aaadara\n")], check_lcs)

P('ABC185_E', ['06_lcs_edit'], r'''
void solv(){
    int N, M; cin >> N >> M;
    vector<ll> A(N), B(M); cin >> A >> B;
    cout << edit_distance(A, B) << '\n';
}''', [("4 3\n1 2 1 3\n1 3 1\n", "2\n"), ("4 6\n1 3 2 4\n1 5 2 6 4 3\n", "3\n"), ("5 5\n1 1 1 1 1\n2 2 2 2 2\n", "5\n")])

P('ABC324_C', ['06_lcs_edit'], r'''
void solv(){
    int N; string T; cin >> N >> T;
    vector<int> ans;
    rep(i, N) { string S; cin >> S; if (within_one_edit(S, T)) ans.push_back(i + 1); }
    cout << sizei(ans) << '\n' << join(ans) << '\n';
}''', [("5 ababc\nababc\nbabc\nabacbc\nabdbc\nabbac\n", "4\n1 2 3 4\n"), ("1 aoki\ntakahashi\n", "0\n\n"),
       ("9 atcoder\natoder\natcode\nathqcoder\natcoder\ntacoder\njttcoder\natoder\natceoder\natcoer\n", "6\n1 2 4 7 8 9\n")])

P('ABC386_C', ['06_lcs_edit'], r'''
void solv(){
    int K; string S, T; cin >> K >> S >> T;
    if (within_one_edit(S, T)) YES; else NO;
}''', [("1\nabc\nagc\n", "Yes\n"), ("1\nabc\nawtf\n", "No\n"), ("1\nabc\nac\n", "Yes\n"), ("1\nback\nblack\n", "Yes\n"), ("1\nsame\nsame\n", "Yes\n"), ("1\nleap\nread\n", "No\n")])

P('ABC386_F', ['06_lcs_edit'], r'''
void solv(){
    int K; string S, T; cin >> K >> S >> T;
    if (edit_distance_bounded(S, T, K) >= 0) YES; else NO;
}''', [("3\nabc\nawtf\n", "Yes\n"), ("2\nabc\nawtf\n", "No\n"), ("17\ntwothousandtwentyfour\nhappynewyear\n", "Yes\n")])

# ================= 06_lis =================
P('ABC038_D', ['06_lis'], r'''
void solv(){
    int N; cin >> N;
    vector<pair<int, int>> b(N); cin >> b;
    cout << sizei(lis_2d(b)) << '\n';
}''', [("3\n3 3\n1 1\n2 2\n", "3\n"), ("2\n4 5\n4 3\n", "1\n"), ("4\n2 5\n3 3\n4 5\n6 6\n", "3\n"), ("5\n8 8\n5 3\n2 2\n4 2\n2 1\n", "4\n")])

P('ABC134_E', ['06_lis'], r'''
void solv(){
    int N; cin >> N;
    vector<ll> A(N); cin >> A;
    cout << min_increasing_cover(A) << '\n';
}''', [("5\n2\n1\n4\n5\n3\n", "2\n"), ("4\n0\n0\n0\n0\n", "4\n")])

P('ABC354_F', ['06_lis'], r'''
void solv(){
    int T; cin >> T;
    while (T--) {
        int N; cin >> N;
        vector<ll> A(N); cin >> A;
        const LisMembers mb = lis_members(A);
        vector<int> idx;
        rep(i, N) if (mb.in_some[i]) idx.push_back(i + 1);
        cout << sizei(idx) << '\n' << join(idx) << '\n';
    }
}''', [("1\n5\n2 1 4 5 3\n", "4\n1 2 3 4\n"), ("2\n6\n2 5 3 4 3 4\n5\n10000 1000 100 1 10\n", "5\n1 3 4 5 6\n2\n4 5\n")])

P('ABC360_G', ['06_lis'], r'''
void solv(){
    int N; cin >> N;
    vector<ll> A(N); cin >> A;
    const auto e = lis_ending_at(A), s = lis_starting_at(A);
    int ans = *max_element(nall(e));
    // 1 か所を変えて、左の列（j）と右の列（k）をつなぐ：k - j ≥ 2、A[k] - A[j] ≥ 2
    rep(k, 1, N) chmax(ans, 1 + s[k]);       // 左は空（k の前を小さい値に）
    rep(j, N - 1) chmax(ans, e[j] + 1);       // 右は空
    vector<ll> vals(A);
    sort(nall(vals)); vals.erase(unique(nall(vals)), vals.end());
    const int K = sizei(vals);
    vector<int> bit(K + 1, 0);  // 値の番号の接頭辞での e[j] の最大
    for (int k = 2; k < N; ++k) {
        const int j = k - 2;
        for (int x = (int)(lower_bound(nall(vals), A[j]) - vals.begin()) + 1; x <= K; x += x & -x) chmax(bit[x], e[j]);
        int best = 0;
        for (int x = (int)(upper_bound(nall(vals), A[k] - 2) - vals.begin()); x > 0; x -= x & -x) chmax(best, bit[x]);
        if (best > 0) chmax(ans, best + 1 + s[k]);
    }
    cout << ans << '\n';
}''', [("4\n3 2 2 4\n", "3\n"), ("5\n4 5 3 6 7\n", "4\n")])

# ================= 07_bracket =================
P('TYP002', ['07_bracket'], r'''
void solv(){
    int N; cin >> N;
    for_each_balanced(N, [](const string& s) { cout << s << '\n'; });
}''', [("2\n", "()\n"), ("3\n", ""), ("4\n", "(())\n()()\n"),
       ("10\n", "((((()))))\n(((()())))\n(((())()))\n(((()))())\n(((())))()\n((()(())))\n((()()()))\n((()())())\n((()()))()\n((())(()))\n((())()())\n((())())()\n((()))(())\n((()))()()\n(()((())))\n(()(()()))\n(()(())())\n(()(()))()\n(()()(()))\n(()()()())\n(()()())()\n(()())(())\n(()())()()\n(())((()))\n(())(()())\n(())(())()\n(())()(())\n(())()()()\n()(((())))\n()((()()))\n()((())())\n()((()))()\n()(()(()))\n()(()()())\n()(()())()\n()(())(())\n()(())()()\n()()((()))\n()()(()())\n()()(())()\n()()()(())\n()()()()()\n")])

P('ABC167_F', ['07_bracket'], r'''
void solv(){
    int N; cin >> N;
    vector<string> S(N); cin >> S;
    if (bracket_concat_order(S).ok) YES; else NO;
}''', [("2\n)\n(()\n", "Yes\n"), ("2\n)(\n()\n", "No\n"), ("4\n((()))\n((((((\n))))))\n()()()\n", "Yes\n"), ("3\n(((\n)\n)\n", "No\n")])

P('ABC223_F', ['07_bracket'], r'''
void solv(){
    int N, Q; string S; cin >> N >> Q >> S;
    BracketSegTree bt(S);
    while (Q--) {
        int t, l, r; cin >> t >> l >> r;
        if (t == 1) {
            const char a = bt.get(l - 1), b = bt.get(r - 1);
            bt.set(l - 1, b), bt.set(r - 1, a);
        } else {
            if (bt.is_balanced(l - 1, r)) YES; else NO;
        }
    }
}''', [("5 3\n(())(\n2 1 4\n2 1 2\n2 4 5\n", "Yes\nNo\nNo\n"), ("5 3\n(())(\n2 1 4\n1 1 4\n2 1 4\n", "Yes\nNo\n"),
       ("8 8\n(()(()))\n2 2 7\n2 2 8\n1 2 5\n2 3 4\n1 3 4\n1 3 5\n1 1 4\n1 6 8\n", "Yes\nNo\nNo\n")])

P('ABC307_D', ['07_bracket'], r'''
void solv(){
    int N; string S; cin >> N >> S;
    cout << erase_matched_pairs(S) << '\n';
}''', [("8\na(b(d))c\n", "ac\n"), ("5\na(b)(\n", "a(\n"), ("2\n()\n", "\n"), ("6\n)))(((\n", ")))(((\n")])

P('ABC350_F', ['07_bracket'], r'''
void solv(){
    string S; cin >> S;
    const auto m = match_brackets(S);
    const auto d = bracket_depths(S);
    string out;
    for (int i = 0, dir = 1; 0 <= i && i < sizei(S); i += dir) {
        if (m[i] >= 0) i = m[i], dir = -dir;  // 組の反対側に飛んで向きを変える
        else out += d[i] % 2 ? swap_case(string(1, S[i]))[0] : S[i];
    }
    cout << out << '\n';
}''', [("((A)y)x\n", "YAx\n"), ("((XYZ)n(X(y)Z))\n", "XYZNXYZ\n"), ("(((()))(()))(())\n", "\n"),
       ("dF(qT(plC())NnnfR(GsdccC))PO()KjsiI((ysA)eWW)ve\n", "dFGsdccCrFNNnplCtQPOKjsiIwwEysAve\n")])

P('ABC394_D', ['07_bracket'], r'''
void solv(){
    string S; cin >> S;
    if (is_balanced(S, "()[]<>")) YES; else NO;
}''', [("([])<>()\n", "Yes\n"), ("([<)]>\n", "No\n"), ("())\n", "No\n")])

P('ABC428_C', ['07_bracket'], r'''
void solv(){
    int Q; cin >> Q;
    vector<BracketBalance> st = {BracketBalance{}};
    while (Q--) {
        int t; cin >> t;
        if (t == 1) { char c; cin >> c; st.push_back(st.back() + balance_of(c)); }
        else st.pop_back();
        if (st.back().balanced()) YES; else NO;
    }
}''', [("8\n1 (\n2\n1 (\n1 )\n2\n1 (\n1 )\n1 )\n", "No\nYes\nNo\nYes\nNo\nNo\nNo\nYes\n")])

# ================= 08_trie =================
P('ABC268_G', ['08_trie', M], r'''
void solv(){
    int N; cin >> N;
    vector<string> S(N); cin >> S;
    Trie TR;
    for (auto& s : S) TR.insert(s);
    using mint = ModInt998;
    const mint half = mint(2).inv();
    for (auto& s : S) {
        const ll P = TR.count_prefixes_of(s) - 1, Q = TR.count_prefix(s) - 1;  // 自分の接頭辞 / 自分で始まる（自分を除く）
        cout << mint(1 + P) + mint(N - 1 - P - Q) * half << '\n';
    }
}''', [("3\na\naa\nab\n", "1\n499122179\n499122179\n"), ("3\na\naa\naaa\n", "1\n2\n3\n")])

P('ABC287_E', ['08_trie'], r'''
void solv(){
    int N; cin >> N;
    vector<string> S(N); cin >> S;
    Trie TR;
    for (auto& s : S) TR.insert(s);
    for (auto& s : S) cout << TR.max_lcp(s, true) << '\n';
}''', [("3\nabc\nabb\naac\n", "2\n2\n1\n"), ("11\nabracadabra\nbracadabra\nracadabra\nacadabra\ncadabra\nadabra\ndabra\nabra\nbra\nra\na\n", "4\n3\n2\n1\n0\n1\n0\n4\n3\n2\n1\n")])

P('ABC353_E', ['08_trie'], r'''
void solv(){
    int N; cin >> N;
    Trie TR;
    rep(i, N) { string s; cin >> s; TR.insert(s); }
    cout << TR.sum_lcp_pairs() << '\n';
}''', [("3\nab abc arc\n", "4\n"), ("11\nab bb aaa bba baba babb aaaba aabbb a a b\n", "32\n")])

P('ABC377_G', ['08_trie'], r'''
void solv(){
    int N; cin >> N;
    Trie TR;
    vector<ll> best;  // best[v] = min(|S_j| - depth(v))
    rep(k, N) {
        string s; cin >> s;
        const ll n = sizei(s);
        ll a = n;
        for (int v : TR.path(s)) if (v < sizei(best)) chmin(a, n - TR.depth(v) + best[v]);
        cout << a << '\n';
        TR.insert(s);
        best.resize(TR.node_count(), LINF);
        for (int v : TR.path(s)) chmin(best[v], n - TR.depth(v));
    }
}''', [("3\nsnuke\nsnuki\nsnuuk\n", "5\n2\n4\n"), ("3\nabc\narc\nagc\n", "3\n3\n3\n"), ("8\nat\natatat\nattat\naatatatt\nattattat\nttatta\ntta\ntt\n", "2\n4\n3\n8\n3\n6\n3\n1\n")])

P('ABC403_E', ['08_trie'], r'''
void solv(){
    int Q; cin >> Q;
    Trie TR;
    vector<char> banned;
    while (Q--) {
        int t; string s; cin >> t >> s;
        if (t == 1) {
            const int v = TR.insert(s, 0);  // ノードだけ作って印を付ける
            banned.resize(TR.node_count(), 0);
            banned[v] = 1;
            TR.erase_prefix(s);
        } else {
            bool ok = true;
            for (int v : TR.path(s)) if (v < sizei(banned) && banned[v]) ok = false;
            if (ok) TR.insert(s);
        }
        cout << TR.size() << '\n';
    }
}''', [("4\n1 at\n2 watcoder\n2 atcoder\n1 wa\n", "0\n1\n1\n0\n"),
       ("10\n1 w\n1 avko\n2 atcoder\n1 bzginn\n2 beginner\n1 atco\n2 contest\n1 ntxcdg\n1 atc\n1 contest\n", "0\n0\n1\n1\n2\n1\n2\n2\n2\n1\n")])

# ================= 09_reversible_deque =================
P('ABC158_D', ['09_reversible_deque'], r'''
void solv(){
    string S; int Q; cin >> S >> Q;
    ReversibleDeque<char> d(S);
    while (Q--) {
        int t; cin >> t;
        if (t == 1) d.reverse();
        else { int f; char c; cin >> f >> c; if (f == 1) d.push_front(c); else d.push_back(c); }
    }
    cout << d.to_string() << '\n';
}''', [("a\n4\n2 1 p\n1\n2 2 c\n1\n", "cpa\n"), ("a\n6\n2 2 a\n2 1 b\n1\n2 2 c\n1\n1\n", "aabc\n"), ("y\n1\n2 1 x\n", "xy\n")])

P('ABC199_C', ['09_reversible_deque'], r'''
void solv(){
    int N, Q; string S; cin >> N >> S >> Q;
    RotatableArray<char> ra(S);
    while (Q--) {
        int t, a, b; cin >> t >> a >> b;
        if (t == 2) ra.rotate_left(N);
        else { const char x = ra[a - 1]; ra[a - 1] = ra[b - 1]; ra[b - 1] = x; }
    }
    cout << ra.to_string() << '\n';
}''', [("2\nFLIP\n2\n2 0 0\n1 1 4\n", "LPFI\n"), ("2\nFLIP\n6\n1 1 3\n2 0 0\n1 1 2\n1 2 3\n2 0 0\n1 1 4\n", "ILPF\n")])

P('ABC237_D', ['09_reversible_deque'], r'''
void solv(){
    int N; string S; cin >> N >> S;
    ReversibleDeque<int> d;
    d.push_back(N);
    for (int i = N; i >= 1; --i) {
        if (S[i - 1] == 'L') d.push_back(i - 1);
        else d.push_front(i - 1);
    }
    cout << join(d.to_vector()) << '\n';
}''', [("5\nLRRLR\n", "1 2 4 5 3 0\n"), ("7\nLLLLLLL\n", "7 6 5 4 3 2 1 0\n")])
