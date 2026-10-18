# 最終照合の広い洗い出し（AtCoder_文字列問題一覧.txt の問題。作ったものを見ていないエージェントが解いたもの）
# agent_2 の担当問題（string_v2 最終照合）
# 実行: cd /home/claude/final/agent_2 && python3 /home/claude/final/run.py probs.py /mnt/user-data/uploads/claude/Main_base.cpp /home/claude/final/out /home/claude/ac-library [問題名…]
MV = '/mnt/user-data/uploads/claude/math_v2/'
GV = '/mnt/user-data/uploads/claude/graph_grid_v2/'

# ============================ C〜F ============================

# ABC282 C String Delimiter：" で区切って偶数番目の , を . に（split ／ replace_all ／ join）
P('ABC282_C', ['02_pattern_match'], r'''
void solv(){
    int n; string s; cin >> n >> s;
    auto v = split(s, '"');                       // 空の部分も残る → 偶数番目が括りの外
    for (size_t i = 0; i < v.size(); i += 2) v[i] = replace_all(v[i], ",", ".");
    cout << join(v, "\"") << '\n';
}''', [
    ('8\n"a,b"c,d\n', '"a,b"c.d\n'),
    ('5\n,,,,,\n', '.....\n'),
    ('20\na,"t,"c,"o,"d,"e,"r,\n', 'a."t,"c."o,"d."e,"r.\n'),
])

# ABC215 C One More aab aba baa：重複のある並べ替えの辞書順 K 番目（math_v2 04_permutation）
P('ABC215_C', [MV + '04_permutation'], r'''
void solv(){
    string s; long long k; cin >> s >> k;
    cout << multiset_perm_unrank(s, k - 1) << '\n';
}''', [
    ('aab 2\n', 'aba\n'),
    ('baba 4\n', 'baab\n'),
    ('ydxwacbz 40320\n', 'zyxwdcba\n'),
])

# ABC381 C 11/22 Substring：ランレングスで '/' 1 個の両隣が 1 の塊・2 の塊か
P('ABC381_C', ['00_str_base'], r'''
void solv(){
    int n; string s; cin >> n >> s;
    auto r = run_length(s);
    long long ans = 1;  // '/' は必ずある
    for (int i = 1; i + 1 < (int)r.size(); i++)
        if (r[i].value == '/' && r[i].len == 1 && r[i - 1].value == '1' && r[i + 1].value == '2')
            ans = max(ans, 1 + 2 * min(r[i - 1].len, r[i + 1].len));
    cout << ans << '\n';
}''', [
    ('8\n211/2212\n', '5\n'),
    ('5\n22/11\n', '1\n'),
    ('22\n/1211/2///2111/2222/11\n', '7\n'),
])

# ABC370 C Word Ladder：下げる位置を左から、上げる位置を右から（std）
P('ABC370_C', ['00_str_base'], r'''
void solv(){
    string s, t; cin >> s >> t;
    int n = s.size();
    vector<string> x;
    rep(i, n) if (s[i] > t[i]) { s[i] = t[i]; x.push_back(s); }
    rrep(i, n) if (s[i] < t[i]) { s[i] = t[i]; x.push_back(s); }
    cout << x.size() << '\n';
    for (auto& e : x) cout << e << '\n';
}''', [
    ('adbe\nbcbc\n', '3\nacbe\nacbc\nbcbc\n'),
    ('abcde\nabcde\n', '0\n'),
    ('afwgebrw\noarbrenq\n', '8\naawgebrw\naargebrw\naarbebrw\naarbebnw\naarbebnq\naarbeenq\naarbrenq\noarbrenq\n'),
])

# ABC407 C Security 2：後ろから B を押す回数 b_i = b_{i+1} + (S_i − S_{i+1}) mod 10（std）
P('ABC407_C', ['00_str_base'], r'''
void solv(){
    string s; cin >> s;
    int n = s.size();
    long long b = s[n - 1] - '0';
    for (int i = n - 2; i >= 0; i--) b += ((s[i] - s[i + 1]) % 10 + 10) % 10;
    cout << n + b << '\n';
}''', [
    ('21\n', '4\n'),
    ('407\n', '17\n'),
    ('2025524202552420255242025524\n', '150\n'),
])

# ABC314 C Rotate Colored Subsequence：色ごとの文字を右に 1 回す（RotatableArray）
P('ABC314_C', ['09_reversible_deque'], r'''
void solv(){
    int n, m; string s; cin >> n >> m >> s;
    vector<vector<int>> pos(m);
    rep(i, n) { int c; cin >> c; pos[c - 1].push_back(i); }
    string t = s;
    for (auto& p : pos) {
        string cs; for (int i : p) cs += s[i];
        RotatableArray<char> ra(cs);
        ra.rotate_right(1);
        rep(j, p.size()) t[p[j]] = ra[j];
    }
    cout << t << '\n';
}''', [
    ('8 3\napzbqrcs\n1 2 3 1 2 2 1 2\n', 'cszapqbr\n'),
    ('2 1\naa\n1 1\n', 'aa\n'),
])

# ABC452 C Fishbones：ok[長さ][位置][文字] を作って各 S_j を確かめる（std）
P('ABC452_C', ['00_str_base'], r'''
void solv(){
    int n; cin >> n;
    vector<int> A(n), B(n);
    rep(i, n) cin >> A[i] >> B[i];
    int m; cin >> m;
    vector<string> S(m); cin >> S;
    static bool ok[11][11][26];
    for (auto& s : S) rep(p, s.size()) ok[s.size()][p][s[p] - 'a'] = true;
    for (auto& s : S) {
        bool g = (int)s.size() == n;
        if (g) rep(i, n) if (!ok[A[i]][B[i] - 1][s[i] - 'a']) { g = false; break; }
        if (g) YES; else NO;
    }
}''', [
    ('5\n5 3\n5 2\n4 1\n5 1\n3 2\n8\nretro\nchris\nitchy\ntuna\ncrab\nrock\ncod\nash\n', 'Yes\nYes\nNo\nNo\nNo\nNo\nNo\nNo\n'),
    ('5\n5 1\n5 2\n5 3\n5 4\n5 5\n8\nretro\nchris\nitchy\ntuna\ncrab\nrock\ncod\nash\n', 'Yes\nYes\nYes\nNo\nNo\nNo\nNo\nNo\n'),
])

# ABC286 C Rotate and Palindrome：回す回数 k ごとに A k + B × min_changes_to_palindrome（O(N^2)）
P('ABC286_C', ['04_palindrome'], r'''
void solv(){
    long long n, a, b; string s; cin >> n >> a >> b >> s;
    long long ans = LINF;
    rep(k, n) chmin(ans, a * k + b * min_changes_to_palindrome(s.substr(k) + s.substr(0, k)));
    cout << ans << '\n';
}''', [
    ('5 1 2\nrrefa\n', '3\n'),
    ('8 1000000000 1000000000\nbcdfcgaa\n', '4000000000\n'),
])

# ABC375 D ABA：真ん中 j ごとに Σ_c（左の c の数）×（右の c の数）（CharPrefix）
P('ABC375_D', ['00_str_base'], r'''
void solv(){
    string s; cin >> s;
    int n = s.size();
    const Alphabet U = Alphabet::upper();
    CharPrefix cp(s, U);
    long long ans = 0;
    rep(j, n) rep(c, 26) ans += (long long)cp.count(0, j, U.ch(c)) * cp.count(j + 1, n, U.ch(c));
    cout << ans << '\n';
}''', [
    ('ABCACC\n', '5\n'),
    ('OOOOOOOO\n', '56\n'),
    ('XYYXYYXYXXX\n', '75\n'),
])

# ABC454 D (xx)：(xx) → xx をスタックで消えなくなるまで（規則が重ならないので結果は一意）→ 一致するか（std）
P('ABC454_D', ['00_str_base'], r'''
string normal_form(const string& s) {
    string st;
    for (char c : s) {
        st.push_back(c);
        size_t k = st.size();
        if (k >= 4 && st.compare(k - 4, 4, "(xx)") == 0) { st.resize(k - 4); st += "xx"; }
    }
    return st;
}
void solv(){
    int T; cin >> T;
    while (T--) {
        string a, b; cin >> a >> b;
        if (normal_form(a) == normal_form(b)) YES; else NO;
    }
}''', [
    ('6\n(xx)x\nx(xx)\n(x)x\n(xx)\n)x()x(\n)x()x(\nx\n(x)\n(((((xx)))))x\nx((((((((((xx))))))))))\n((xx)xx)xx\n(x((xx))x)(xx)\n',
     'Yes\nNo\nYes\nNo\nYes\nYes\n'),
])

# ABC418 D XNOR Operation：美しい ⇔ 0 の数が偶数 → 先頭からの 0 の数の偶奇が同じ組を数える（std）
P('ABC418_D', ['00_str_base'], r'''
void solv(){
    int n; string t; cin >> n >> t;
    long long cnt[2] = {1, 0}, ans = 0;
    int p = 0;
    for (char c : t) { if (c == '0') p ^= 1; ans += cnt[p]++; }
    cout << ans << '\n';
}''', [
    ('3\n110\n', '3\n'),
    ('4\n0000\n', '4\n'),
    ('30\n011011100101110111100010011010\n', '225\n'),
])

# ABC426 D Pop and Insert：残す塊（run）ごとに N + （その文字の総数）− 2 × 長さ の最小（run_length ／ count_chars）
P('ABC426_D', ['00_str_base'], r'''
void solv(){
    int T; cin >> T;
    while (T--) {
        int n; string s; cin >> n >> s;
        auto cnt = count_chars(s, Alphabet::digits());
        long long ans = LINF;
        for (auto& r : run_length(s)) chmin(ans, n + cnt[r.value - '0'] - 2 * r.len);
        cout << ans << '\n';
    }
}''', [
    ('3\n5\n01001\n3\n000\n15\n110010111100101\n', '4\n0\n16\n'),
])

# ABC312 D Count Bracket Sequences：開きの数 − 閉じの数 を状態にする DP（math_v2 ModInt998）
P('ABC312_D', [MV + '00_modint'], r'''
void solv(){
    string s; cin >> s;
    int n = s.size();
    using mint = ModInt998;
    vector<mint> dp(n + 2, 0);
    dp[0] = 1;
    for (char c : s) {
        vector<mint> nd(n + 2, 0);
        rep(j, n + 1) {
            if (c != ')') nd[j + 1] += dp[j];
            if (c != '(' && j > 0) nd[j - 1] += dp[j];
        }
        dp = nd;
    }
    cout << dp[0] << '\n';
}''', [
    ('(???(?\n', '2\n'),
    (')))))\n', '0\n'),
    ('??????????????(????????(??????)?????????(?(??)\n', '603032273\n'),
])

# ABC455 E Unbalanced ABC Substrings：全体 −（a=b・b=c・c=a の和 − 2 × 全部等しい）（先頭からの差を map で数える。std）
P('ABC455_E', ['00_str_base'], r'''
void solv(){
    int n; string s; cin >> n >> s;
    map<long long, long long> m1, m2, m3;
    map<pair<int, int>, long long> m4;
    int a = 0, b = 0, c = 0;
    long long e = 0, e4 = 0;
    auto add = [&] { e += m1[a - b]++ + m2[b - c]++ + m3[c - a]++; e4 += m4[{a - b, b - c}]++; };
    add();
    for (char ch : s) { (ch == 'A' ? a : ch == 'B' ? b : c)++; add(); }
    cout << (long long)n * (n + 1) / 2 - (e - 2 * e4) << '\n';
}''', [
    ('6\nAABBCC\n', '4\n'),
    ('6\nABCABC\n', '0\n'),
    ('10\nACABCAABAB\n', '17\n'),
])

# ABC473 F：A を開き・B を閉じとして、区間に「対応しない閉じ」がないか（BracketSegTree(S, "AB")）
P('ABC473_F', ['07_bracket'], r'''
void solv(){
    int n; string s; cin >> n >> s;
    BracketSegTree bt(s, "AB");
    int q; cin >> q;
    while (q--) {
        int t; cin >> t;
        if (t == 1) { int i; char c; cin >> i >> c; bt.set(i - 1, c); }
        else { int l, r; cin >> l >> r; if (bt.balance(l - 1, r).need_open == 0) YES; else NO; }
    }
}''', [
    ('10\nAABBAABABB\n6\n2 1 10\n1 5 B\n2 1 10\n2 6 8\n1 3 A\n2 1 10\n', 'Yes\nNo\nYes\nYes\n'),
])

# ABC359 D Avoid K Palindrome：直前 K−1 文字のビットを状態にする DP。長さ K の回文の表は is_palindrome で
P('ABC359_D', [MV + '00_modint', '00_str_base'], r'''
void solv(){
    int n, k; string s; cin >> n >> k >> s;
    using mint = ModInt998;
    const int W = 1 << (k - 1);
    vector<char> pal(1 << k);
    rep(m, 1 << k) {
        string t; rep(j, k) t += (m >> j & 1) ? 'B' : 'A';
        pal[m] = is_palindrome(t);
    }
    auto fits = [&](int i, int x) { return s[i] == '?' || s[i] == (x ? 'B' : 'A'); };
    vector<mint> dp(W, 0);
    rep(m, W) {  // 最初の K−1 文字（新しい文字ほど下のビット）
        bool ok = true;
        rep(j, k - 1) ok = ok && fits(j, m >> (k - 2 - j) & 1);
        if (ok) dp[m] = 1;
    }
    rep(i, k - 1, n) {
        vector<mint> nd(W, 0);
        rep(m, W) rep(x, 2) {
            if (!fits(i, x)) continue;
            const int full = (int)(m << 1 | x);
            if (!pal[full]) nd[full & (W - 1)] += dp[m];
        }
        dp = nd;
    }
    mint ans = 0;
    for (auto& v : dp) ans += v;
    cout << ans << '\n';
}''', [
    ('7 4\nAB?A?BA\n', '1\n'),
    ('40 7\n????????????????????????????????????????\n', '116295436\n'),
    ('15 5\nABABA??????????\n', '0\n'),
    ('40 8\n?A?B??B?B?AA?A?B??B?A???B?BB?B???BA??BAA\n', '259240\n'),
])

# ============================ B ============================

# ABC202 B 180°：反転して 6 と 9 を入れ替える（reversed ／ CharMap）
P('ABC202_B', ['00_str_base'], r'''
void solv(){
    string s; cin >> s;
    CharMap cm(Alphabet::digits());
    cm.set('6', '9'); cm.set('9', '6');
    cout << cm.apply(reversed(s)) << '\n';
}''', [
    ('0601889\n', '6881090\n'),
    ('86910\n', '01698\n'),
    ('01010\n', '01010\n'),
])

# ABC290 B Qual B：K 個目より後の o を x に（std）
P('ABC290_B', ['00_str_base'], r'''
void solv(){
    int n, k; string s; cin >> n >> k >> s;
    for (char& c : s) if (c == 'o') { if (k > 0) k--; else c = 'x'; }
    cout << s << '\n';
}''', [
    ('10 3\noxxoxooxox\n', 'oxxoxoxxxx\n'),
])

# ABC426 B The Odd One Out：1 個だけの文字（count_chars）
P('ABC426_B', ['00_str_base'], r'''
void solv(){
    string s; cin >> s;
    auto cnt = count_chars(s);
    rep(c, 26) if (cnt[c] == 1) cout << Alphabet::lower().ch(c) << '\n';
}''', [
    ('odd\n', 'o\n'),
    ('dad\n', 'a\n'),
    ('wwwwwwwwwv\n', 'v\n'),
])

# ABC380 B Hurdle Parsing：| で区切って長さ（split ／ join）
P('ABC380_B', ['00_str_base'], r'''
void solv(){
    string s; cin >> s;
    vector<int> a;
    for (auto& p : split(s, '|', true)) a.push_back(p.size());
    cout << join(a) << '\n';
}''', [
    ('|---|-|----|-|-----|\n', '3 1 4 1 5\n'),
    ('|----------|\n', '10\n'),
    ('|-|-|-|------|\n', '1 1 1 6\n'),
])

# ABC394 B cat：長さの順に並べてつなげる（sort ＋ join(v, "")）
P('ABC394_B', ['00_str_base'], r'''
void solv(){
    int n; cin >> n;
    vector<string> v(n); cin >> v;
    sort(nall(v), [](const string& a, const string& b) { return a.size() < b.size(); });
    cout << join(v, "") << '\n';
}''', [
    ('3\ntc\noder\na\n', 'atcoder\n'),
    ('4\ncat\nenate\non\nc\n', 'concatenate\n'),
])

# ABC459 B：先頭の文字を電話のキーの数字に（std の表）
P('ABC459_B', ['00_str_base'], r'''
void solv(){
    const string key = "22233344455566677778889999";
    int n; cin >> n;
    string ans;
    rep(i, n) { string w; cin >> w; ans += key[w[0] - 'a']; }
    cout << ans << '\n';
}''', [
    ('2\nalgorithm heuristic\n', '24\n'),
    ('3\ni love you\n', '459\n'),
])

# ABC279 B LOOKUP：T が S の部分文字列か（find_first）
P('ABC279_B', ['02_pattern_match'], r'''
void solv(){
    string s, t; cin >> s >> t;
    if (find_first(s, t) != -1) YES; else NO;
}''', [
    ('voltage\ntag\n', 'Yes\n'),
    ('atcoder\nace\n', 'No\n'),
    ('gorilla\ngorillagorillagorilla\n', 'No\n'),
    ('toyotasystems\ntoyotasystems\n', 'Yes\n'),
])

# ABC367 B Cut .0：末尾の 0 と . を消す（std）
P('ABC367_B', ['00_str_base'], r'''
void solv(){
    string x; cin >> x;
    while (x.back() == '0') x.pop_back();
    if (x.back() == '.') x.pop_back();
    cout << x << '\n';
}''', [
    ('1.012\n', '1.012\n'),
    ('12.340\n', '12.34\n'),
    ('99.900\n', '99.9\n'),
    ('0.000\n', '0\n'),
])

# ABC287 B Postal Card：下 3 桁が T のどれかか（std set ＋ substr）
P('ABC287_B', ['00_str_base'], r'''
void solv(){
    int n, m; cin >> n >> m;
    vector<string> s(n), t(m); cin >> s >> t;
    set<string> st(nall(t));
    int ans = 0;
    for (auto& x : s) ans += st.count(x.substr(3));
    cout << ans << '\n';
}''', [
    ('3 3\n142857\n004159\n071028\n159\n287\n857\n', '2\n'),
    ('5 4\n235983\n109467\n823476\n592801\n000333\n333\n108\n467\n983\n', '3\n'),
    ('4 4\n000000\n123456\n987111\n000000\n000\n111\n999\n111\n', '3\n'),
])

# ABC338 B Frequency：最も多い文字（count_chars ＋ max_element。同数なら番号の小さい方）
P('ABC338_B', ['00_str_base'], r'''
void solv(){
    string s; cin >> s;
    auto cnt = count_chars(s);
    cout << Alphabet::lower().ch(max_element(nall(cnt)) - cnt.begin()) << '\n';
}''', [
    ('frequency\n', 'e\n'),
    ('atcoder\n', 'a\n'),
    ('pseudopseudohypoparathyroidism\n', 'o\n'),
])

# ABC413 B cat 2：S_i + S_j（i ≠ j）の種類（std set）
P('ABC413_B', ['00_str_base'], r'''
void solv(){
    int n; cin >> n;
    vector<string> s(n); cin >> s;
    set<string> st;
    rep(i, n) rep(j, n) if (i != j) st.insert(s[i] + s[j]);
    cout << st.size() << '\n';
}''', [
    ('4\nat\natco\ncoder\nder\n', '11\n'),
    ('5\na\naa\naaa\naaaa\naaaaa\n', '7\n'),
    ('10\narmiearggc\nukupaunpiy\ncogzmjmiob\nrtwbvmtruq\nqapfzsitbl\nvhkihnipny\nybonzypnsn\nesxvgoudra\nusngxmaqpt\nyfseonwhgp\n', '90\n'),
])

# ABC354 B AtCoder Janken 2：名前を辞書順に並べて T mod N 番目（std）
P('ABC354_B', ['00_str_base'], r'''
void solv(){
    int n; cin >> n;
    vector<string> name(n);
    long long t = 0;
    rep(i, n) { long long c; cin >> name[i] >> c; t += c; }
    sort(nall(name));
    cout << name[t % n] << '\n';
}''', [
    ('3\ntakahashi 2\naoki 6\nsnuke 5\n', 'snuke\n'),
    ('3\ntakahashi 2813\ntakahashixx 1086\ntakahashix 4229\n', 'takahashix\n'),
])

# ABC421 B Fibonacci Reversed：和を文字列にして反転（reversed）
P('ABC421_B', ['00_str_base'], r'''
void solv(){
    long long x, y; cin >> x >> y;
    vector<long long> a = {x, y};
    rep(i, 8) a.push_back(stoll(reversed(to_string(a[i] + a[i + 1]))));
    cout << a[9] << '\n';
}''', [
    ('1 1\n', '415\n'),
    ('3 7\n', '895\n'),
    ('90701 90204\n', '9560800101\n'),
])

# ABC277 B Playing Cards Validation：文字の種類の判定（Alphabet::of(…).contains）＋ 重複（set）
P('ABC277_B', ['00_str_base'], r'''
void solv(){
    int n; cin >> n;
    vector<string> s(n); cin >> s;
    const Alphabet suit = Alphabet::of("HDCS"), rank = Alphabet::of("A23456789TJQK");
    bool ok = set<string>(nall(s)).size() == (size_t)n;
    for (auto& x : s) ok = ok && suit.contains(x[0]) && rank.contains(x[1]);
    if (ok) YES; else NO;
}''', [
    ('4\nH3\nDA\nD3\nSK\n', 'Yes\n'),
    ('5\nH3\nDA\nCK\nH3\nS7\n', 'No\n'),
    ('4\n3H\nAD\n3D\nKS\n', 'No\n'),
    ('5\n00\nAA\nXX\nYY\nZZ\n', 'No\n'),
])

# ABC281 B Sandwich Number：大文字・100000〜999999・大文字（std）
P('ABC281_B', ['00_str_base'], r'''
void solv(){
    string s; cin >> s;
    bool ok = s.size() == 8 && isupper(s[0]) && isupper(s[7]) && s[1] != '0';
    rep(i, 1, 7) ok = ok && isdigit(s[i]);
    if (ok) YES; else NO;
}''', [
    ('Q142857Z\n', 'Yes\n'),
    ('AB912278C\n', 'No\n'),
    ('X900000\n', 'No\n'),
    ('K012345K\n', 'No\n'),
])

# ABC403 B Four Hidden：? つきの T に U が現れるか（wildcard_match）
P('ABC403_B', ['02_wildcard_match'], r'''
void solv(){
    string t, u; cin >> t >> u;
    if (!wildcard_match(t, u).empty()) YES; else NO;
}''', [
    ('tak??a?h?\nnashi\n', 'Yes\n'),
    ('??e??e\nsnuke\n', 'No\n'),
    ('????\naoki\n', 'Yes\n'),
])

# ABC349 B Commencement：出現回数ごとの文字の種類が 0 か 2 か（count_chars ＋ map）
P('ABC349_B', ['00_str_base'], r'''
void solv(){
    string s; cin >> s;
    map<int, int> f;
    for (int c : count_chars(s)) if (c > 0) f[c]++;
    bool ok = true;
    for (auto& [k, v] : f) ok = ok && v == 2;
    if (ok) YES; else NO;
}''', [
    ('commencement\n', 'Yes\n'),
    ('banana\n', 'No\n'),
    ('ab\n', 'Yes\n'),
])

# ABC285 B Longest Uncommon Prefix：ずらし i ごとに先頭から違う文字が続く長さ（std。O(N^2)）
P('ABC285_B', ['00_str_base'], r'''
void solv(){
    int n; string s; cin >> n >> s;
    rep(i, 1, n) {
        int l = 0;
        while (l + i < n && s[l] != s[l + i]) l++;
        cout << l << '\n';
    }
}''', [
    ('6\nabcbac\n', '5\n1\n2\n0\n1\n'),
])

# ABC366 B Vertical Writing：* で長さをそろえて時計回りに 90 度（graph_grid_v2 rotate_cw）→ 末尾の * を消す
P('ABC366_B', ['00_str_base', GV + '10_grid_transform'], r'''
void solv(){
    int n; cin >> n;
    vector<string> s(n); cin >> s;
    size_t m = 0;
    for (auto& x : s) m = max(m, x.size());
    for (auto& x : s) x.resize(m, '*');
    for (auto r : rotate_cw(s)) {
        while (r.back() == '*') r.pop_back();
        cout << r << '\n';
    }
}''', [
    ('3\nabc\nde\nfghi\n', 'fda\ngeb\nh*c\ni\n'),
    ('3\natcoder\nbeginner\ncontest\n', 'cba\noet\nngc\ntio\nend\nsne\nter\n*r\n'),
])
