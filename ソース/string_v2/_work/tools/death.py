#!/usr/bin/env python3
"""使えない引数でメッセージを出して止まるかを確かめる（作業用）
使い方: death.py <文字列の out フォルダ>"""
import os, subprocess, sys
out = sys.argv[1]
START = 'ここからテンプレート本体（コピー範囲）'
END = '// ======================= ここまでテンプレート本体'
def body(name):
    t = open(os.path.join(out, name + '.cpp'), encoding='utf-8-sig').read()
    head = t[:t.index(START)]
    inc = '\n'.join(l for l in head.split('\n') if l.startswith('#include'))
    return inc + '\n' + t[t.index('\n', t.index(START)) + 1:t.index(END)]
CASES = [
    # (ファイル, 文, stderr に含まれるべき文字列)
    ('01_rolling_hash', 'RollingHash H("abc"); H.get(2, 1);', 'RollingHash'),
    ('01_rolling_hash', 'RollingHash H("abc"); H.get(0, 4);', 'RollingHash'),
    ('02_pattern_match', 'find_first("abc", "b", 5);', 'find_first'),
    ('03_suffix_array', 'SuffixArray S("abc"); S.kth_substring(100);', 'kth_substring'),
    ('03_suffix_array', 'SuffixArray S("abc"); S.lcp(0, 5);', 'SuffixArray'),
    ('04_palindrome', 'Palindromes P("abc"); P.is_palindrome(1, 0);', 'Palindromes'),
    ('05_lex_order', 'dict_unrank(-1, 3);', 'dict_unrank'),
    ('06_subsequence', 'NextTable nx("abc"); nx.match_length("a", 5);', 'NextTable'),
    ('06_subsequence', 'NextTable nx("aBc");', 'Alphabet::id'),
    ('06_subsequence', 'count_distinct_subsequences<long long>(std::string(80, \'a\') + std::string("bcdefghijklmnopqrstuvwxyzbcdefghijklmnopqrstuvwxyzbcdefghijklmnopqrstuvwxyz"));', 'count_distinct_subsequences'),
    ('06_subsequence', 'count_subsequence_occurrences<long long>(std::string(67, \'a\') + "b", std::string(33, \'a\') + "b");', 'count_subsequence_occurrences'),
    ('04_palindrome', 'count_palindromic_subsequences<long long>(std::string(64, \'a\'));', 'count_palindromic_subsequences'),
    ('06_lcs_edit', 'edit_distance_bounded("abc", "abd", -1);', 'edit_distance_bounded'),
    ('06_lcs_edit', 'edit_distance("ab", "b", -1, 1, 1);', 'edit_distance'),
    ('06_lcs_edit', 'shortest_common_supersequence(std::string("a"), std::vector<int>{1000});', 'shortest_common_supersequence'),
    ('06_lis', 'std::vector<long long> v; for (int i = 0; i < 64; ++i) v.push_back(2 * (i / 1) + 0), v.push_back(2 * i - 1); count_lis<long long>(v);', 'count_lis'),
    ('07_bracket', 'BracketSegTree bt(std::vector<int>{296, 41});', 'BracketSegTree'),
    ('08_trie', 'Trie T; T.insert("a", 4500000000LL); T.sum_lcp_pairs();', 'Trie::sum_lcp_pairs'),
    ('03_suffix_array', 'auto S = SuffixArray::of_strings(std::vector<std::string>{"ab", "b"}); S.offset(2);', 'SuffixArray::offset'),
    ('03_suffix_array', 'SuffixArray S("ab"); S.value(2);', 'SuffixArray::value'),
    ('03_suffix_array', 'SuffixArray S("ab"); S.position(1, 0);', 'SuffixArray::position'),
    ('06_lcs_edit', 'edit_script(std::string(20000, \'a\'), std::string(20000, \'b\'));', 'edit_script'),
    ('06_lis', 'std::vector<int> a; for (int i = 0; i < 200; ++i) a.push_back(i / 2 * 2 + (i % 2 ? 0 : 1)); count_lis<long long>(a);', 'count_lis'),
    ('07_bracket', 'is_balanced("()", "(");', 'is_balanced'),
    ('07_bracket', 'is_balanced("()", "(())");', 'is_balanced'),
    ('07_bracket', 'balance_of("()", "()[]");', 'balance_of'),
    ('07_bracket', 'BracketSegTree bt("()[]", "()[]");', 'BracketSegTree'),
    ('07_bracket', 'BracketSegTree bt("()"); bt.set(2, \'(\');', 'BracketSegTree::set'),
    ('07_bracket', 'StaticBracket sb("()"); sb.is_balanced(0, 3);', 'StaticBracket'),
    ('08_trie', 'Trie T; T.insert("abC");', 'Alphabet::id'),
    ('08_trie', 'Trie T; T.insert("ab"); T.kth(1);', 'Trie::kth'),
    ('08_trie', 'Trie T; T.insert("ab"); T.max_lcp("a", true);', 'Trie::max_lcp'),
    ('08_trie', 'Trie T; T.insert("ab", -1);', 'Trie::insert'),
    ('08_trie', 'Trie T; T.parent(5);', 'Trie::parent'),
    ('09_reversible_deque', 'ReversibleDeque<char> d; d.pop_back();', 'ReversibleDeque::pop_back'),
    ('09_reversible_deque', 'ReversibleDeque<char> d("ab"); d[2];', 'ReversibleDeque::operator[]'),
    ('09_reversible_deque', 'ReversibleDeque<char> d("ab"); d.substr(1, 2);', 'ReversibleDeque::substr'),
    ('09_reversible_deque', 'RotatableArray<int> ra(std::vector<int>{1, 2}); ra.set(2, 0);', 'RotatableArray::set'),
    # Phase 4
    ('01_set_hash', 'MultisetPrefix mp(std::vector<int>{1, 2}); mp.get(1, 3);', 'MultisetPrefix::get'),
    ('02_aho_corasick', 'AhoCorasick ac; ac.add("ab"); ac.build(); ac.add("b");', 'AhoCorasick::add'),
    ('02_aho_corasick', 'AhoCorasick ac; ac.add("");', 'AhoCorasick::add'),
    ('02_aho_corasick', 'AhoCorasick ac; ac.add("ab"); ac.count_total("abab");', 'AhoCorasick'),
    ('02_aho_corasick', 'AhoCorasick ac; ac.add("ab"); ac.build(); ac.next(5, \'a\');', 'AhoCorasick::next'),
    ('02_aho_corasick', 'AhoCorasick ac; ac.add("ab"); ac.build(); ac.node_of(1);', 'AhoCorasick::node_of'),
    ('02_aho_corasick', 'AhoCorasick ac; ac.add("aB");', 'Alphabet::id'),
    ('02_wildcard_match', 'wildcard_match(std::string(9000000, \'a\'), std::string(1, \'a\'));', 'wildcard_match'),
    ('02_wildcard_match', 'mismatch_counts(std::vector<int>(9000000), std::vector<int>(1));', 'mismatch_counts'),
    ('03_suffix_automaton', 'SuffixAutomaton S("abc"); S.kth_substring(6);', 'SuffixAutomaton::kth_substring'),
    ('03_suffix_automaton', 'SuffixAutomaton S("abc"); S.kth_substring(-1, true);', 'SuffixAutomaton::kth_substring'),
    ('03_suffix_automaton', 'SuffixAutomaton S("abc"); S.len(100);', 'SuffixAutomaton::len'),
    ('03_suffix_automaton', 'SuffixAutomaton S("abc"); S.extend(\'Z\');', 'Alphabet::id'),
    ('04_palindromic_tree', 'PalindromicTree P("abc"); P.link(10);', 'PalindromicTree::link'),
    ('04_palindromic_tree', 'PalindromicTree P("abc"); P.node_of_prefix(4);', 'PalindromicTree::node_of_prefix'),
    ('04_palindromic_tree', 'PalindromicTree P(std::string(100, \'a\')); P.count_partitions<long long>();', 'PalindromicTree::count_partitions'),
    ('09_mutable_string', 'MutableString M("abc"); M.set(3, \'a\');', 'MutableString::set'),
    ('09_mutable_string', 'MutableString M("abc"); M.kth_char(0, 2, 2);', 'MutableString::kth_char'),
    ('09_mutable_string', 'MutableString M("abc"); M.hash(2, 1);', 'MutableString::hash'),
    ('09_mutable_string', 'MutableString M("abc"); M.lcp(0, 4);', 'MutableString::lcp'),
    ('09_mutable_string', 'MutableString M(std::vector<int>{300}, Alphabet::ints(301)); M.to_string();', 'MutableString::to_string'),
    ('03_suffix_automaton', 'SuffixAutomaton S("abc"); S.next(0, 0);', 'Alphabet::id'),
    ('03_suffix_automaton', 'SuffixAutomaton S(Alphabet::ints(200000000)); S.extend(1);', 'SuffixAutomaton'),
    ('02_aho_corasick', 'AhoCorasick ac; ac.add("ab"); ac.build(); ac.contained_patterns(-1);', 'AhoCorasick::contained_patterns'),
    ('02_aho_corasick', 'AhoCorasick ac(Alphabet::ints(300000000));', 'AhoCorasick'),
    ('04_palindromic_tree', 'PalindromicTree P(Alphabet::ints(1000000000));', 'PalindromicTree'),
    ('04_palindromic_tree', 'PalindromicTree P(std::string(100, \'a\')); P.count_partitions_prefix<long long>();', 'PalindromicTree::count_partitions_prefix'),
    ('01_set_hash', 'MultisetPrefix mp(std::vector<int>{1, 2}); mp.same(0, 3, 0, 1);', 'MultisetPrefix::same'),
    ('09_mutable_string', 'MutableString M("abc"); M.same(0, 100, 0, 1);', 'MutableString::same'),
    ('09_mutable_string', 'MutableString M(std::vector<int>(1000000), Alphabet::ints(300)); M.count(0, 1, 0);', 'MutableString'),
    ('00_str_base', 'CharPrefix cp(std::string(1000000, \'a\'), Alphabet::bytes());', 'CharPrefix'),
    # 最終照合
    ('02_pattern_match', 'rewrite_until_stable(std::string("ab"), "a", "aa");', 'rewrite_until_stable'),
    ('02_pattern_match', 'rewrite_until_stable(std::string("ab"), "", "x");', 'rewrite_until_stable'),
    ('02_pattern_match', 'Occurrences oc("abc", "b"); oc.count(0, 4);', 'Occurrences::count'),
    ('00_str_base', 'from_ids(std::vector<long long>{26});', 'from_ids'),
    ('00_str_base', 'CharPrefix cp("abc"); cp.count_id(0, 3, 26);', 'CharPrefix::count_id'),
    ('00_str_base', 'split("a,b", std::string(""));', 'split'),
    ('00_str_base', 'Alphabet::ints(3).with("a");', 'Alphabet::with'),
    ('00_str_base', 'Alphabet::lower().with("a");', 'Alphabet::of'),
    ('04_palindrome', 'count_palindromes_leq<long long>(std::string(30, \'z\'));', 'count_palindromes_leq'),
    ('05_lex_order', 'lex_rank_same_length<long long>(std::string(14, \'z\'));', 'lex_rank_same_length'),
    ('07_bracket', 'count_balanced_completions<long long>(std::string("(a)"));', 'count_balanced_completions'),
    ('07_bracket', 'count_balanced_completions<long long>(std::string(130, \'?\'));', 'count_balanced_completions'),
    ('01_rolling_hash', 'center_palindrome_lengths_with_changes(std::string("ab"), -1);', 'center_palindrome_lengths_with_changes'),
]
ok = True
cache = {}
for i, (name, stmt, want) in enumerate(CASES):
    if name not in cache: cache[name] = body(name)
    src = cache[name] + '\nint main() { ' + stmt + ' std::cout << "NOT STOPPED\\n"; return 0; }\n'
    fn = f'/tmp/death_{i}.cpp'
    open(fn, 'w', encoding='utf-8').write(src)
    r = subprocess.run(['g++', '-std=gnu++17', '-O1', '-w', '-o', f'/tmp/death_{i}', fn], capture_output=True, text=True)
    if r.returncode:
        print('COMPILE FAIL', name, stmt, r.stderr[:800]); ok = False; continue
    p = subprocess.run([f'/tmp/death_{i}'], capture_output=True, text=True)
    good = p.returncode != 0 and 'NOT STOPPED' not in p.stdout and want in p.stderr
    ok &= good
    print('OK ' if good else 'NG ', name, '|', stmt[:70], '|', p.stderr.strip().split('\n')[0][:120])
print('ALL OK' if ok else 'FAILED')
