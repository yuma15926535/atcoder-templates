void solv(){
    // 文字列第二世代（基盤）
    string S = "abracadabra";
    CharPrefix cp(S);
    CharPositions cps(S);
    auto runs = run_length(S);
    auto ids = to_ids(S);
    auto cv = compress_values(vector<ll>{5, -3, 5});
    debug(cp.count(0, 4, 'a'), cps.next(1, 'a'), runs, StrSpan{1, 3}, ids, cv.ids);
    // ACL の文字列の関数も今までどおり呼べる
    auto sa = suffix_array(S);
    auto lcp = lcp_array(S, sa);
    auto z = z_algorithm(S);
    cout << join(sa) << "\n" << join(lcp) << "\n" << join(z) << "\n";
    // 数学の is_palindrome_str と文字列の is_palindrome が共存
    cout << is_palindrome_str("level") << is_palindrome("level") << "\n";
    cout << reversed(S) << " " << min_with_reverse(S) << " " << split("a,b", ',').size() << "\n";
    vector<int> v = {3, 1, 2};
    cout << sorted_chars(v)[0] << "\n";
    Alphabet A = Alphabet::upper();
    cout << A.id('C') << "\n";
    CharMap cm;
    cm.replace('a', 'b');
    cout << cm.apply(S) << "\n";
    using mint = ModInt998;
    vector<mint> mv = {mint(3), mint(-1)};
    cout << join(mv) << "\n";
    // グラフ・グリッドの trim（グリッド）と一緒に
    vector<string> g = {"..", ".#"};
    auto tg = trim(g, '.');
    cout << tg.size() << " " << count_chars(S)[0] << " " << hamming_distance("ab", "ba") << " " << common_prefix_length(S, "abc") << "\n";
    // 文字列第二世代（Phase 3）
    {
        RollingHash H(S);
        cout << H.lcp(0, 7) << " " << H.same(0, 4, 7, 11) << " " << (H.get(0, 1) == hash_of("a")) << "\n";
        auto zz = z_array(S);
        auto pi = prefix_function(S);
        cout << find_all(S, "abra").size() << " " << count_occurrences(S, "a") << " " << zz[7] << " " << pi.back() << " " << min_period(S) << "\n";
        SuffixArray SAx(S);
        auto sa2 = make_suffix_array(S);
        auto lcp2 = make_lcp_array(S, sa2);
        cout << SAx.count_distinct_substrings() << " " << (sa2 == sa) << " " << (lcp2 == lcp) << " " << SAx.count("abra") << "\n";
        Palindromes P(S);
        cout << P.longest().len() << " " << P.count_all() << " " << min_changes_to_palindrome(S) << "\n";
        cout << min_rotation(S) << " " << concat_order(vector<string>{"b", "ba"}).size() << " " << dict_rank("ab", 3) << "\n";
        NextTable nx(S);
        cout << nx.is_subsequence("aaa") << " " << count_distinct_subsequences<mint>(S) << " " << shortest_non_subsequence(S) << "\n";
        cout << lcs_length(S, "cadabra") << " " << edit_distance("kitten", "sitting") << " " << lcs_restore(S, "abc") << "\n";
        vector<ll> a = {3, 1, 4, 1, 5, 9, 2, 6};
        auto lis = lis_length(a);
        cout << lis << " " << join(lis_restore(a)) << " " << count_lis<mint>(a) << " " << min_increasing_cover(a) << "\n";
        cout << is_balanced("(()())") << " " << join(match_brackets("(())")) << " " << balance_of(")(") << "\n";
        BracketSegTree bt("(())");
        bt.flip(0);
        cout << bt.all_balanced() << " " << StaticBracket("()()").is_balanced(1, 3) << "\n";
        Trie trie;
        for (string w : {"ab", "abc", "b"}) trie.insert(w);
        cout << trie.count_prefix("ab") << " " << trie.kth(1) << " " << trie.rank("b") << " " << trie.sum_lcp_pairs() << "\n";
        ReversibleDeque<char> dq("abc");
        dq.reverse();
        dq.push_front('x');
        RotatableArray<int> ra(vector<int>{1, 2, 3});
        ra.rotate_left(1);
        cout << dq.to_string() << " " << join(ra.to_vector()) << "\n";
        int rank = 0, hash = 0, next = 0, prev = 0, count = 0, runs = 0, periods = 0, lcs = 0, sa = 0, lcp = 0, rev = 0;  // よく使う変数名
        cout << rank + hash + next + prev + count + runs + periods + lcs + sa + lcp + rev << "\n";
    }
    // 文字列第二世代（Phase 4）
    {
        vector<int> a = {1, 2, 1, 3};
        auto ph = prefix_set_hashes(a);
        MultisetPrefix mp(a);
        SetHash sh;
        sh.add(5);
        MultisetHash mh;
        mh.add(5, 2);
        cout << (ph[2] == ph[3]) << " " << mp.same(0, 2, 1, 3) << " " << sh.distinct() << " " << mh.count << " " << (value_hash(3) != 0) << "\n";
        AhoCorasick ac;
        for (string w : {"he", "she", "his", "hers"}) ac.add(w);
        ac.build();
        cout << ac.count_total("ushers") << " " << join(ac.count_each("ushers")) << " " << ac.contains_any("xyz") << "\n";
        cout << join(wildcard_match("abcab", "a?")) << " " << join(mismatch_counts("abcab", "ab")) << " " << match_with_k_mismatches("abcab", "cb", 1).size() << "\n";
        SuffixAutomaton sam(S);
        cout << sam.count_distinct_substrings() << " " << sam.count_occurrences("abra") << " " << sam.kth_substring(3) << " " << sam.longest_common_substring("cadabra").len() << "\n";
        PalindromicTree pt(S);
        cout << pt.count_distinct() << " " << pt.count_total() << " " << pt.min_partition() << " " << pt.count_partitions<mint>() << "\n";
        cout << lyndon_factorization(S).size() << " " << is_lyndon("aab") << " " << join(lyndon_array("abaab")) << " " << maximal_repetitions(S).size() << " "
             << count_square_substrings("aaaa") << "\n";
        MutableString ms(S);
        ms.set(0, 'z');
        cout << ms.to_string() << " " << ms.lcp(0, 7) << " " << ms.is_palindrome(0, 1) << " " << ms.kth_char(0, 4, 0) << " " << ms.count(0, 11, 'a') << "\n";
        int link = 0, len = 0, occ = 0, cnt = 0, fail = 0, go = 0, wild = 0, lyndon = 0, word = 0, period = 0, node = 0, state = 0;  // よく使う変数名
        cout << link + len + occ + cnt + fail + go + wild + lyndon + word + period + node + state << "\n";
    }
    // 文字列第二世代（最終照合で足したもの）
    {
        cout << rewrite_until_stable(string("WWA"), "WA", "AC") << " " << erase_until_stable(string("ABCABC"), "ABC").size() << " " << Occurrences(S, "abra").count(0, 11) << "\n";
        cout << rearrange_no_adjacent_equal(string("aab")) << " " << lex_smallest_no_adjacent_equal(string("aabb")) << " " << split_runs_by(S, [](char a, char b) { return a != b; }).size() << "\n";
        cout << subsequence_positions("ace", "abcde").size() << " " << common_prefix_length("a?c", "abd", '?') << " " << strip_right("ab**", '*') << " " << join(split("a, b", ", "), '|') << "\n";
        cout << count_palindromes_leq<mint>("AXA", Alphabet::upper()) << " " << lex_rank_same_length<ll>("ba") << " " << count_balanced_completions<mint>(string("(??)")) << "\n";
        cout << center_palindrome_lengths_with_changes(S, 1).size() << " " << to_upper('a') << caesar_shift('z', 1) << " " << from_ids(vector<ll>{0, 1}) << " " << count_chars(S)[0] << "\n";
        int strip = 0, rewrite = 0, occ = 0, runs2 = 0;  // よく使いそうな変数名
        cout << strip + rewrite + occ + runs2 << "\n";
    }
    rep(i, 3) cout << caesar_shift("a", i);
    cout << "\nsolv done\n";
}
