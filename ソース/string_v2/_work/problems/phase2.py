# Phase 2（00_str_base）の代表的な問題
P('ABC259_C', ['00_str_base'], r'''
void solv(){
    string S, T; cin >> S >> T;
    auto a = run_length(S), b = run_length(T);
    bool ok = a.size() == b.size();
    for (size_t i = 0; ok && i < a.size(); ++i)
        ok = a[i].value == b[i].value && (a[i].len == b[i].len || (a[i].len >= 2 && a[i].len <= b[i].len));
    if (ok) YES; else NO;
}''', [("abbaac\nabbbbaaac\n", "Yes\n"), ("xyzz\nxyyzz\n", "No\n")])

P('ABC329_C', ['00_str_base'], r'''
void solv(){
    int N; string S; cin >> N >> S;
    vector<ll> best(26, 0);
    for (auto& r : run_length(S)) chmax(best[r.value - 'a'], r.len);
    cout << accumulate(nall(best), 0LL) << '\n';
}''', [("6\naaabaa\n", "4\n"), ("1\nx\n", "1\n"), ("12\nssskkyskkkky\n", "8\n")])

P('ABC342_C', ['00_str_base'], r'''
void solv(){
    int N, Q; string S; cin >> N >> S >> Q;
    CharMap cm;
    while (Q--) { char c, d; cin >> c >> d; cm.replace(c, d); }
    cout << cm.apply(S) << '\n';
}''', [("7\natcoder\n4\nr a\nt e\nd v\na r\n", "recover\n"), ("3\nabc\n4\na a\ns k\nn n\nz b\n", "abc\n"),
       ("34\nsupercalifragilisticexpialidocious\n20\ng c\nl g\ng m\nc m\nr o\ns e\na a\no f\nf s\ne t\nt l\nd v\np k\nv h\nx i\nh n\nn j\ni r\ns i\nu a\n", "laklimamriiamrmrllrmlrkramrjimrial\n")])

P('ABC310_C', ['00_str_base'], r'''
void solv(){
    int N; cin >> N;
    set<string> st;
    rep(i, N) { string s; cin >> s; st.insert(min_with_reverse(s)); }
    cout << st.size() << '\n';
}''', [("6\na\nabc\nde\ncba\nde\nabc\n", "3\n")])

P('TYP084', ['00_str_base'], r'''
void solv(){
    ll N; string S; cin >> N >> S;
    ll ans = N * (N + 1) / 2;
    for (auto& r : run_length(S)) ans -= r.len * (r.len + 1) / 2;
    cout << ans << '\n';
}''', [("4\nooxo\n", "5\n"), ("5\noxoxo\n", "10\n"), ("5\nooooo\n", "0\n"), ("7\nxxoooxx\n", "16\n")])

P('ABC380_D', ['00_str_base'], r'''
void solv(){
    string S; int Q; cin >> S >> Q;
    const ll n = sizei(S);
    vector<char> ans;
    rep(i, Q) {
        ll K; cin >> K; --K;
        const char c = S[K % n];
        ans.push_back(popcnt(K / n) & 1 ? swap_case(string(1, c))[0] : c);
    }
    cout << join(ans) << '\n';
}''', [("aB\n16\n1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16\n", "a B A b A b a B A b a B a B A b\n"),
       ("qWeRtYuIoP\n8\n1 1 2 3 5 8 13 21\n", "q q W e t I E Q\n"),
       ("AnUoHrjhgfLMcDIpzxXmEWPwBZvbKqQuiJTtFSlkNGVReOYCdsay\n5\n1000000000000000000 123456789 1 987654321 999999999999999999\n", "K a A Z L\n")])

P('ABC232_B', ['00_str_base'], r'''
void solv(){
    string S, T; cin >> S >> T;
    bool ok = false;
    rep(k, 26) ok |= caesar_shift(S, k) == T;
    if (ok) YES; else NO;
}''', [("abc\nijk\n", "Yes\n"), ("z\na\n", "Yes\n"), ("ppq\nqqp\n", "No\n"), ("atcoder\natcoder\n", "Yes\n")])

P('ABC219_C', ['00_str_base'], r'''
void solv(){
    string X; int N; cin >> X >> N;
    const Alphabet A = Alphabet::of(X);   // X の順が辞書順
    vector<string> S(N);
    cin >> S;
    sort(nall(S), [&](const string& a, const string& b) { return to_ids(a, A) < to_ids(b, A); });
    for (auto& s : S) cout << s << '\n';
}''', [("bacdefghijklmnopqrstuvwxzy\n4\nabx\nbzz\nbzy\ncaa\n", "bzz\nbzy\nabx\ncaa\n"),
       ("zyxwvutsrqponmlkjihgfedcba\n5\na\nab\nabc\nac\nb\n", "b\na\nac\nab\nabc\n")])
