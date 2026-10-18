# Phase 4（☆ 7 本）の代表的な問題。ABC の入出力例は AtCoder の問題ページで確かめたもの
# Library Checker（LC_）は入力例が公式（library-checker-problems の gen/example_*.in）、出力は下の愚直で作ったもの

# ---------------- Library Checker の出力を愚直で作る ----------------
def _is_lyndon(w): return len(w) > 0 and all(w < w[k:] for k in range(1, len(w)))
def bf_lyndon_factorization(s):
    cut, i = [0], 0
    while i < len(s):
        best = max(l for l in range(1, len(s) - i + 1) if _is_lyndon(s[i:i + l]))  # いちばん長い Lyndon 語の接頭辞
        i += best
        cut.append(i)
    return ' '.join(map(str, cut)) + '\n'
def _period(w):
    return next(p for p in range(1, len(w) + 1) if all(w[i] == w[i + p] for i in range(len(w) - p)))
def bf_runs(s):
    n, res = len(s), []
    for l in range(n):
        for r in range(l + 2, n + 1):
            p = _period(s[l:r])
            if 2 * p > r - l: continue
            if l > 0 and s[l - 1] == s[l - 1 + p]: continue
            if r < n and s[r] == s[r - p]: continue
            res.append((p, l, r))
    res.sort()
    return f'{len(res)}\n' + ''.join(f'{p} {l} {r}\n' for p, l, r in res)
def bf_eertree(s):
    pal = lambda w: w == w[::-1]
    num, order = {}, []
    for j in range(1, len(s) + 1):
        for i in range(j):
            w = s[i:j]
            if pal(w) and w not in num: num[w] = len(order) + 1; order.append(w)
    lines = [str(len(order))]
    for w in order:
        p = -1 if len(w) == 1 else (0 if len(w) == 2 else num[w[1:-1]])
        suf = [w[k:] for k in range(1, len(w)) if pal(w[k:])]
        lines.append(f'{p} {num[suf[0]] if suf else 0}')
    v = []
    for i in range(1, len(s) + 1):
        w = next(s[k:i] for k in range(i) if pal(s[k:i]))
        v.append(str(num[w]))
    lines.append(' '.join(v))
    return '\n'.join(lines) + '\n'
def bf_distinct(s): return f'{len({s[l:r] for l in range(len(s)) for r in range(l + 1, len(s) + 1)})}\n'

# ================= 01_set_hash =================
P('ABC250_E', ['01_set_hash'], r'''
void solv(){
    int N; cin >> N;
    vector<ll> a(N), b(N);
    rep(i, N) cin >> a[i];
    rep(i, N) cin >> b[i];
    const auto ha = prefix_set_hashes(a), hb = prefix_set_hashes(b);  // 先頭 i 個の集合
    int Q; cin >> Q;
    while (Q--) {
        int x, y; cin >> x >> y;
        if (ha[x] == hb[y]) YES; else NO;
    }
}''', [("5\n1 2 3 4 5\n1 2 2 4 3\n7\n1 1\n2 2\n2 3\n3 3\n4 4\n4 5\n5 5\n", "Yes\nYes\nYes\nNo\nNo\nYes\nNo\n")])

P('ABC367_F', ['01_set_hash'], r'''
void solv(){
    int N, Q; cin >> N >> Q;
    vector<int> a(N), b(N);
    rep(i, N) cin >> a[i];
    rep(i, N) cin >> b[i];
    const MultisetPrefix A(a), B(b);  // 値の乱数は A と B で共通なので、別の列どうしでも比べられる
    while (Q--) {
        int l, r, L, R; cin >> l >> r >> L >> R;
        if (r - l == R - L && A.get(l - 1, r) == B.get(L - 1, R)) YES; else NO;
    }
}''', [("5 4\n1 2 3 2 4\n2 3 1 4 2\n1 3 1 3\n1 2 3 5\n1 4 2 5\n1 5 1 5\n", "Yes\nNo\nNo\nYes\n"),
       ("4 4\n4 4 4 4\n4 4 4 4\n1 2 2 3\n3 3 1 1\n1 3 1 4\n1 4 2 3\n", "Yes\nYes\nNo\nNo\n")])

# ================= 02_aho_corasick =================
P('ABC268_Ex', ['02_aho_corasick'], r'''
void solv(){
    string S; int N; cin >> S >> N;
    AhoCorasick ac;
    rep(i, N) { string t; cin >> t; ac.add(t); }
    ac.build();
    // 左から読み、どれかのパターンが終わった瞬間にその文字を * にする（状態を根に戻す）のが最適
    int v = 0, ans = 0;
    for (char c : S) {
        v = ac.next(v, c);
        if (ac.match_count(v) > 0) ++ans, v = 0;
    }
    cout << ans << '\n';
}''', [("abcdefghijklmn\n3\nabcd\nijk\nghi\n", "2\n"), ("atcoderbeginnercontest\n1\nabc\n", "0\n"), ("aaaaaaaaa\n2\naa\nxyz\n", "4\n")])

P('ABC354_G', ['02_aho_corasick'], r'''
void solv(){
    int N; cin >> N;
    vector<string> S(N);
    vector<ll> A(N);
    rep(i, N) cin >> S[i];
    rep(i, N) cin >> A[i];
    map<string, ll> best;  // 同じ文字列はどれか 1 つしか選べないので、重みが最大のものだけ残す
    rep(i, N) best[S[i]] = max(best[S[i]], A[i]);
    AhoCorasick ac;
    vector<ll> w;
    for (auto& [s, a] : best) ac.add(s), w.push_back(a);
    ac.build();
    // 含む関係（推移的）で重み最大の反鎖 = 重みの合計 − 最大流（i が j を含むなら 左 i → 右 j）
    const int n = (int)w.size(), src = 2 * n, snk = 2 * n + 1;
    mf_graph<ll> g(2 * n + 2);
    ll total = 0;
    rep(i, n) {
        g.add_edge(src, i, w[i]), g.add_edge(n + i, snk, w[i]), total += w[i];
        for (int j : ac.contained_patterns(i)) g.add_edge(i, n + j, LINF);
    }
    cout << total - g.flow(src, snk) << '\n';
}''', [("4\natcoder\nat\ncoder\ncode\n5 2 3 4\n", "6\n"), ("10\nabcd\nabc\nab\na\nb\nc\nd\nab\nbc\ncd\n100 10 50 30 60 90 80 70 40 20\n", "260\n")])

P('ABC362_G_AC', ['02_aho_corasick'], r'''
void solv(){
    string S; int Q; cin >> S >> Q;
    AhoCorasick ac;
    rep(i, Q) { string t; cin >> t; ac.add(t); }
    ac.build();
    for (ll c : ac.count_each(S)) cout << c << '\n';
}''', [("missisippi\n5\ni\ns\na\nis\nmissisippi\n", "4\n3\n0\n2\n1\n"), ("aaaaaa\n6\na\naa\naaa\naaaa\naaaaa\naaaaaa\n", "6\n5\n4\n3\n2\n1\n")])

# ================= 02_wildcard_match =================
P('ABC196_F', ['02_wildcard_match'], r'''
void solv(){
    string S, T; cin >> S >> T;
    const auto mc = mismatch_counts(S, T);  // ずらしごとの不一致の数（σ = 2 なので畳み込み 2 回分）
    cout << *min_element(mc.begin(), mc.end()) << '\n';
}''', [("0001\n101\n", "1\n"), ("0101010\n1010101\n", "7\n"), ("10101000010011011110\n0010011111\n", "1\n")])

# ================= 03_suffix_automaton =================
P('ABC362_G_SAM', ['03_suffix_automaton'], r'''
void solv(){
    string S; int Q; cin >> S >> Q;
    const SuffixAutomaton sam(S);
    while (Q--) {
        string t; cin >> t;
        cout << sam.count_occurrences(t) << '\n';
    }
}''', [("missisippi\n5\ni\ns\na\nis\nmissisippi\n", "4\n3\n0\n2\n1\n"), ("aaaaaa\n6\na\naa\naaa\naaaa\naaaaa\naaaaaa\n", "6\n5\n4\n3\n2\n1\n")])

_ns = ['abcbcba', 'mississippi']
P('LC_number_of_substrings', ['03_suffix_automaton'], r'''
void solv(){
    string S; cin >> S;
    cout << SuffixAutomaton(S).count_distinct_substrings() << '\n';
}''', [(s + '\n', bf_distinct(s)) for s in _ns])

# ================= 04_palindromic_tree =================
P('ABC237_Ex', ['04_palindromic_tree'], r'''
void solv(){
    string S; cin >> S;
    const PalindromicTree pt(S);
    const int m = pt.size(), n = m - 2;  // 回文はノード 2〜m-1
    // v に含まれる回文 = v から link（回文の接尾辞）と parent（両端を除く）をたどって着くもの
    // 含まれない最大の選び方（反鎖）= 回文の数 − 最大マッチング
    const int src = 2 * n, snk = 2 * n + 1;
    mf_graph<int> g(2 * n + 2);
    for (int v = 2; v < m; ++v) {
        vector<char> seen(m, 0);
        vector<int> st = {v};
        while (!st.empty()) {
            const int u = st.back();
            st.pop_back();
            for (int w : {pt.link(u), pt.parent(u)})
                if (w >= 2 && !seen[w]) seen[w] = 1, st.push_back(w), g.add_edge(v - 2, n + w - 2, 1);
        }
        g.add_edge(src, v - 2, 1), g.add_edge(n + v - 2, snk, 1);
    }
    cout << n - g.flow(src, snk) << '\n';
}''', [("ababb\n", "3\n"), ("xyz\n", "3\n"), ("xxxxxxxxxx\n", "1\n")])

_es = ['abaa', 'aaaaaaa']
P('LC_eertree', ['04_palindromic_tree'], r'''
void solv(){
    string S; cin >> S;
    const PalindromicTree pt(S);
    // 番号：ODD = -1（ノード 0）、EVEN = 0（ノード 1）、回文は最初に現れた順に 1〜（ノード v は v - 1）
    cout << pt.count_distinct() << '\n';
    for (int v = 2; v < pt.size(); ++v) cout << (pt.len(v) == 1 ? -1 : pt.parent(v) - 1) << ' ' << pt.link(v) - 1 << '\n';
    rep(i, 1, (ll)S.size() + 1) cout << pt.node_of_prefix((int)i) - 1 << " \n"[i == (ll)S.size()];
}''', [(s + '\n', bf_eertree(s)) for s in _es])

# ================= 05_lyndon =================
_ls = ['babaabaab', 'ababacaca']
P('LC_lyndon_factorization', ['05_lyndon'], r'''
void solv(){
    string S; cin >> S;
    cout << 0;
    for (const StrSpan& f : lyndon_factorization(S)) cout << ' ' << f.r;
    cout << '\n';
}''', [(s + '\n', bf_lyndon_factorization(s)) for s in _ls])

_rs = ['abcbcba', 'mississippi', 'ababacaca', 'aaaaa']
P('LC_runenumerate', ['05_lyndon'], r'''
void solv(){
    string S; cin >> S;
    auto runs = maximal_repetitions(S);  // (l, r) の順なので、問題の (周期, l, r) の順に並べ直す
    sort(runs.begin(), runs.end(), [](const Repetition& a, const Repetition& b) { return tie(a.period, a.l, a.r) < tie(b.period, b.l, b.r); });
    cout << runs.size() << '\n';
    for (const Repetition& x : runs) cout << x.period << ' ' << x.l << ' ' << x.r << '\n';
}''', [(s + '\n', bf_runs(s)) for s in _rs])

# ================= 09_mutable_string =================
P('ABC285_F', ['09_mutable_string'], r'''
void solv(){
    int N, Q; string S; cin >> N >> S >> Q;
    MutableString ms(S);
    while (Q--) {
        int t; cin >> t;
        if (t == 1) {
            int x; char c; cin >> x >> c;
            ms.set(x - 1, c);
        } else {
            int l, r; cin >> l >> r; --l;
            // 並べ替えた全体に含まれる ⇔ 区間が昇順で、区間の最小と最大の間の文字は全体にある分を全部含む
            bool ok = ms.is_sorted(l, r);
            if (ok) {
                const auto in = ms.counts(l, r), all = ms.counts(0, N);
                const int lo = ms.alphabet().id(ms.get(l)), hi = ms.alphabet().id(ms.get(r - 1));
                for (int c = lo + 1; c < hi; ++c) ok &= in[c] == all[c];
            }
            if (ok) YES; else NO;
        }
    }
}''', [("6\nabcdcf\n4\n2 1 3\n2 2 6\n1 5 e\n2 2 6\n", "Yes\nNo\nYes\n")])

P('ABC331_F', ['09_mutable_string'], r'''
void solv(){
    int N, Q; string S; cin >> N >> Q >> S;
    MutableString ms(S);
    while (Q--) {
        int t; cin >> t;
        if (t == 1) {
            int x; char c; cin >> x >> c;
            ms.set(x - 1, c);
        } else {
            int L, R; cin >> L >> R;
            if (ms.is_palindrome(L - 1, R)) YES; else NO;
        }
    }
}''', [("7 8\nabcbacb\n2 1 5\n2 4 7\n2 2 2\n1 5 c\n2 1 5\n2 4 7\n1 4 c\n2 3 6\n", "Yes\nNo\nYes\nNo\nYes\nYes\n")])
