# 最終照合（典型 5-5 の残り・洗い出しで足した機能）の問題。入出力例は AtCoder の問題ページで確かめたもの
M = '/mnt/user-data/uploads/claude/math_v2/'
G = '/mnt/user-data/uploads/claude/graph_grid_v2/'

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

# ================= 00_str_base =================
P('ABC242_D', ['00_str_base'], r'''
void solv(){
    string S; int Q; cin >> S >> Q;
    while (Q--) {
        ll t, k; cin >> t >> k; --k;
        ll sh = 0;  // 親（k / 2 文字目）をたどる。子は親を (k % 2 + 1) ずらした文字
        while (t > 0 && k > 0) sh += k % 2 + 1, k /= 2, --t;
        cout << caesar_shift(S[k], sh + t, Alphabet::of("ABC")) << '\n';
    }
}''', [("ABC\n4\n0 1\n1 1\n1 3\n1 6\n", "A\nB\nC\nB\n"),
       ("CBBAACCCCC\n5\n57530144230160008 659279164847814847\n29622990657296329 861239705300265164\n509705228051901259 994708708957785197\n176678501072691541 655134104344481648\n827291290937314275 407121144297426665\n", "A\nA\nC\nA\nA\n")])

P('ABC459_D', ['00_str_base'], r'''
void solv(){
    int T; cin >> T;
    while (T--) {
        string S; cin >> S;
        const string r = rearrange_no_adjacent_equal(S);
        if (r.empty()) { NO; continue; }
        YES;
        cout << r << '\n';
    }
}''', [("3\naiiw\ndoodoo\naabbababcacababaaba\n", "Yes\niwai\nNo\nYes\nababacabababacababa\n")], chk_459d)

P('ABC287_D', ['00_str_base'], r'''
void solv(){
    string S, T; cin >> S >> T;
    const int m = T.size();
    const int p = common_prefix_length(S, T, '?');                             // 前から ? 込みで一致する長さ
    const int q = common_suffix_length(S, T, '?');                             // 後ろから
    rep(x, m + 1) { if (x <= p && m - x <= q) YES; else NO; }
}''', [("a?c\nb?\n", "Yes\nNo\nNo\n"), ("atcoder\n?????\n", "Yes\nYes\nYes\nYes\nYes\nYes\n"),
       ("beginner\ncontest\n", "No\nNo\nNo\nNo\nNo\nNo\nNo\nNo\n")])

P('ABC352_B', ['00_str_base'], r'''
void solv(){
    string S, T; cin >> S >> T;
    auto pos = subsequence_positions(S, T);
    for (auto& p : pos) ++p;
    cout << join(pos) << '\n';
}''', [("abc\naxbxyc\n", "1 3 6\n"), ("aaaa\nbbbbaaaa\n", "5 6 7 8\n"), ("atcoder\natcoder\n", "1 2 3 4 5 6 7\n")])

P('ABC456_C', ['00_str_base'], r'''
void solv(){
    string S; cin >> S;
    ll ans = 0;
    for (const StrSpan& sp : split_runs_by(S, [](char a, char b) { return a != b; })) ans += (ll)sp.len() * (sp.len() + 1) / 2;
    cout << ans % 998244353 << '\n';
}''', [("abbc\n", "6\n"),
       ("cabcabcbcaccacbcbcaabacbacaabccacbccbcacbacbacabcacabcaccaaaaabababcbabacaccabbcacbcbcbcababcbcbabca\n", "760\n")])

P('ABC367_B', ['00_str_base'], r'''
void solv(){
    string X; cin >> X;
    cout << strip_right(strip_right(X, '0'), '.') << '\n';  // 末尾の 0 を除き、残った末尾の . も除く
}''', [("1.012\n", "1.012\n"), ("12.340\n", "12.34\n"), ("99.900\n", "99.9\n"), ("0.000\n", "0\n")])

# ================= 01_rolling_hash =================
P('ABC468_D', ['01_rolling_hash'], r'''
void solv(){
    string S; cin >> S;
    const auto L = center_palindrome_lengths_with_changes(S, 1);  // 1 文字まで書き換えてよい
    ll ans = 0;
    rep(c, L.size()) ans += c % 2 == 0 ? (L[c] + 1) / 2 : L[c] / 2;  // 中心ごとに、長さ L 以下で同じ偶奇のもの
    cout << ans << '\n';
}''', [("ababa\n", "13\n"), ("atcoder\n", "18\n"), ("abccbacbacb\n", "40\n")])

# ================= 02_pattern_match =================
P('ABC394_C', ['02_pattern_match'], r'''
void solv(){
    string S; cin >> S;
    cout << rewrite_until_stable(S, "WA", "AC") << '\n';
}''', [("WACWA\n", "ACCAC\n"), ("WWA\n", "ACC\n"), ("WWWWW\n", "WWWWW\n")])

P('ABC328_D', ['02_pattern_match'], r'''
void solv(){
    string S; cin >> S;
    cout << erase_until_stable(S, "ABC") << '\n';
}''', [("BAABCBCCABCAC\n", "BCAC\n"), ("ABCABC\n", "\n"), ("AAABCABCABCAABCABCBBBAABCBCCCAAABCBCBCC\n", "AAABBBCCC\n")])

P('ABC454_D', ['02_pattern_match'], r'''
void solv(){
    int T; cin >> T;
    while (T--) {
        string A, B; cin >> A >> B;  // (xx) → xx を繰り返した形（順によらない）が同じか
        if (rewrite_until_stable(A, "(xx)", "xx") == rewrite_until_stable(B, "(xx)", "xx")) YES; else NO;
    }
}''', [('6\n(xx)x\nx(xx)\n(x)x\n(xx)\n)x()x(\n)x()x(\nx\n(x)\n(((((xx)))))x\nx((((((((((xx))))))))))\n((xx)xx)xx\n(x((xx))x)(xx)\n',
        'Yes\nNo\nYes\nNo\nYes\nYes\n')])

P('ABC477_C', ['02_pattern_match'], r'''
void solv(){
    int Q; string S, T; cin >> Q >> S >> T;
    const Occurrences occ(S, T);
    while (Q--) {
        int L, R; cin >> L >> R;
        if (occ.exists(L - 1, R)) YES; else NO;
    }
}''', [("3\nabcdabc\nbc\n2 6\n3 5\n3 7\n", "Yes\nNo\nYes\n"),
       ("10\nbbbaaabbbaaaaaabbbab\nbbb\n5 16\n5 6\n2 16\n6 7\n17 19\n3 7\n7 16\n7 20\n17 19\n2 6\n",
        "Yes\nNo\nYes\nNo\nNo\nNo\nYes\nYes\nNo\nNo\n")])

# ================= 03_suffix_array（5-5 の残り） =================
P('ABC240_Ex', ['03_suffix_array'], r'''
void solv(){
    int N; string S; cin >> N >> S;
    // i 番目の部分文字列は長さ i 以下にできる（短くしても増える列のまま）→ 長さは B（B(B+1)/2 > N）以下だけ
    int B = 1;
    while ((ll)B * (B + 1) / 2 <= N) ++B;
    const SuffixArray SA(S);
    // 部分文字列 (l, len) の辞書順のキー = (LCP が len 以上で続く SA の区間の左端, len)。同じ文字列は同じキー
    vector<int> left(N);           // len ごとに更新：left[k] = 順位 k の接尾辞の区間の左端
    vector<ull> keys;
    keys.reserve((size_t)N * B);
    for (int len = 1; len <= B; ++len) {
        for (int k = 0; k < N; ++k) left[k] = (k > 0 && SA.lcp_adj[k - 1] >= len) ? left[k - 1] : k;
        for (int k = 0; k < N; ++k) {
            const int l = SA.sa[k];
            if (l + len <= N) keys.push_back(((ull)left[k] << 40) | ((ull)len << 20) | (ull)l);
        }
    }
    sort(keys.begin(), keys.end());
    vector<int> bit(N + 1, 0);   // 終わりの位置 r での最大の K（先頭からの最大を持つフェニック木）
    auto upd = [&](int i, int v) { for (++i; i <= N; i += i & -i) bit[i] = max(bit[i], v); };
    auto qry = [&](int r) { int v = 0; for (; r > 0; r -= r & -r) v = max(v, bit[r]); return v; };  // [0, r)
    int ans = 0;
    vector<pair<int, int>> pend;
    for (size_t a = 0; a < keys.size();) {
        size_t b = a;
        while (b < keys.size() && (keys[b] >> 20) == (keys[a] >> 20)) ++b;  // 同じ文字列の出現をまとめて
        pend.clear();
        for (size_t t = a; t < b; ++t) {
            const int l = keys[t] & ((1 << 20) - 1), len = (keys[t] >> 20) & ((1 << 20) - 1);
            const int v = qry(l) + 1;  // l より前に終わる、より小さい文字列の続き
            pend.push_back({l + len - 1, v});
            ans = max(ans, v);
        }
        for (auto [r, v] : pend) upd(r, v);
        a = b;
    }
    cout << ans << '\n';
}''', [("7\n0101010\n", "3\n"), ("30\n000011001110101001011110001001\n", "9\n")])

# ================= 04_palindrome（5-5 の残り・足した機能） =================
P('ABC242_E', ['04_palindrome', M + '00_modint'], r'''
using mint = ModInt998;
void solv(){
    int T; cin >> T;
    while (T--) {
        int N; string S; cin >> N >> S;
        cout << count_palindromes_leq<mint>(S, Alphabet::upper()) << '\n';
    }
}''', [("5\n3\nAXA\n6\nABCZAZ\n30\nQWERTYUIOPASDFGHJKLZXCVBNMQWER\n28\nJVIISNEOXHSNEAAENSHXOENSIIVJ\n31\nKVOHEEMSOZZASHENDIGOJRTJVMVSDWW\n",
        "24\n29\n212370247\n36523399\n231364016\n")])

P('ABC349_G', ['04_palindrome', G + '03_union_find'], r'''
void solv(){
    int N; cin >> N;
    vector<int> A(N);
    rep(i, N) cin >> A[i];
    UnionFind uf(N);
    bool ok = true;
    // Manacher と同じ順に、まだ分かっていない組だけ等しくする（合わせて O(N) 回）
    int c = 0, R = -1;
    for (int i = 0; i < N; ++i) {
        int k = 0;
        if (i <= R) {
            const int j = 2 * c - i;
            k = min(A[j], R - i);
            if (A[i] < k || (k < R - i && A[i] != k)) ok = false;  // 対称な位置と矛盾
        }
        for (int d = k + 1; d <= A[i]; ++d) uf.unite(i - d, i + d);
        if (i + A[i] > R) c = i, R = i + A[i];
    }
    vector<vector<int>> diff(N);  // 1 つ外は違う
    rep(i, N) {
        const int l = i - A[i] - 1, r = i + A[i] + 1;
        if (l >= 0 && r < N) {
            if (uf.same(l, r)) ok = false;
            diff[uf.find(l)].push_back(uf.find(r)), diff[uf.find(r)].push_back(uf.find(l));
        }
    }
    vector<int> val(N, 0), S(N);
    rep(i, N) {
        const int r = uf.find(i);
        if (val[r] == 0) {  // まだ決まっていない組には、違うべき組が使っていない最小の数
            set<int> used;
            for (int x : diff[r]) if (val[x]) used.insert(val[x]);
            int v = 1;
            while (used.count(v)) ++v;
            val[r] = v;
        }
        S[i] = val[r];
    }
    // 作った S の各中心の最長の回文の半径が A と同じか（04_palindrome の Manacher で確かめる）
    if (ok) {
        const auto L = center_palindrome_lengths(S);
        rep(i, N) if ((L[2 * i] - 1) / 2 != A[i]) ok = false;
    }
    if (!ok) { NO; return; }
    YES;
    cout << join(S) << '\n';
}''', [("7\n0 0 2 0 2 0 0\n", "Yes\n1 1 2 1 1 1 2\n"), ("7\n0 1 2 3 2 1 0\n", "Yes\n1 1 1 1 1 1 1\n"), ("7\n0 1 2 0 2 1 0\n", "No\n")])

# ================= 07_bracket =================
P('ABC312_D', ['07_bracket', M + '00_modint'], r'''
using mint = ModInt998;
void solv(){
    string S; cin >> S;
    cout << count_balanced_completions<mint>(S) << '\n';
}''', [("(???(?\n", "2\n"), (")))))\n", "0\n"), ("??????????????(????????(??????)?????????(?(??)\n", "603032273\n")])

P('ABC473_F', ['07_bracket'], r'''
void solv(){
    int n; string s; cin >> n >> s;
    BracketSegTree bt(s, "AB");  // A を開き、B を閉じとみる
    int q; cin >> q;
    while (q--) {
        int t; cin >> t;
        if (t == 1) { int i; char c; cin >> i >> c; bt.set(i - 1, c); }
        else { int l, r; cin >> l >> r; if (bt.balance(l - 1, r).need_open == 0) YES; else NO; }  // 途中で B が多くならない
    }
}''', [('10\nAABBAABABB\n6\n2 1 10\n1 5 B\n2 1 10\n2 6 8\n1 3 A\n2 1 10\n', 'Yes\nNo\nYes\nYes\n')])

# ================= 03_suffix_array（5-5：ABC362 G の 3 つ目の解き方） =================
P('ABC362_G_SA', ['03_suffix_array'], r'''
void solv(){
    string S; int Q; cin >> S >> Q;
    const SuffixArray SA(S);
    while (Q--) {
        string t; cin >> t;
        cout << SA.count(t) << '\n';  // O(|t| log N)
    }
}''', [("missisippi\n5\ni\ns\na\nis\nmissisippi\n", "4\n3\n0\n2\n1\n"), ("aaaaaa\n6\na\naa\naaa\naaaa\naaaaa\naaaaaa\n", "6\n5\n4\n3\n2\n1\n")])
