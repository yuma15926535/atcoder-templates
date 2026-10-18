#!/usr/bin/env python3
"""文字列第二世代テンプレートの .cpp を組み立てる（作業用）。
src/<name>.head.txt   … 先頭のコメントと #include
src/<name>.body.hpp   … テンプレート本体（#ifndef GG_XXX ... #endif // GG_XXX）
src/<name>.rest.cpp   … テスト・速度計測・デモ・main（{COMMON_TEST} {COMMON_DEMO} {SETUP_CONSOLE} を差し込む）
依存する本体を順に同梱して out/<番号>_<name>.cpp を作る。
"""
import os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, 'src')
OUT = os.path.join(HERE, 'out')

# 名前: (番号, 依存（同梱する順）, テスト・デモの接頭辞)
FILES = {
    'str_base': ('00', [], 'strb'),
    # Phase 3（★）
    'rolling_hash': ('01', ['str_base'], 'rh'),
    'pattern_match': ('02', ['str_base'], 'pm'),
    'suffix_array': ('03', ['str_base'], 'sa'),
    'palindrome': ('04', ['str_base'], 'pal'),
    'lex_order': ('05', ['str_base'], 'lex'),
    'subsequence': ('06', ['str_base'], 'subq'),
    'lcs_edit': ('06', ['str_base'], 'lcs'),
    'lis': ('06', ['str_base'], 'lis'),
    'bracket': ('07', ['str_base'], 'brk'),
    'trie': ('08', ['str_base'], 'trie'),
    'reversible_deque': ('09', ['str_base'], 'rdq'),
    # Phase 4（☆）
    'set_hash': ('01', ['str_base'], 'shash'),
    'aho_corasick': ('02', ['str_base'], 'ac'),
    'wildcard_match': ('02', ['str_base'], 'wc'),
    'suffix_automaton': ('03', ['str_base'], 'sam'),
    'palindromic_tree': ('04', ['str_base'], 'pt'),
    'lyndon': ('05', ['str_base', 'suffix_array'], 'lyn'),
    'mutable_string': ('09', ['str_base', 'rolling_hash'], 'ms'),
}

START = '// ======================= ここからテンプレート本体（コピー範囲） ======================='
END = '// ======================= ここまでテンプレート本体 ======================='

SETUP_CONSOLE = '''// =============================================================================
//  Windows のコンソールで日本語が文字化けしないようにする（Windows 以外では何もしない）
// =============================================================================
#ifdef _WIN32
extern "C" __declspec(dllimport) int __stdcall SetConsoleOutputCP(unsigned int);
#endif
static void setup_console() {
#ifdef _WIN32
    const char* probe = "あ";  // UTF-8 なら 0xE3 0x81 0x82、CP932 なら 0x82 0xA0
    SetConsoleOutputCP((unsigned char)probe[0] == 0xE3 ? 65001u : 932u);
#endif
}'''


def read(p):
    with open(p, encoding='utf-8') as f:
        return f.read()


def body(name):
    b = read(os.path.join(SRC, name + '.body.hpp')).rstrip('\n')
    guard = 'GG_' + name.upper()
    assert b.startswith('#ifndef ' + guard + '\n#define ' + guard + '\n'), name
    assert b.endswith('#endif  // ' + guard), name
    return b


def build(name):
    num, deps, pre = FILES[name]
    head = read(os.path.join(SRC, name + '.head.txt')).rstrip('\n')
    rest = read(os.path.join(SRC, name + '.rest.cpp'))
    ct = read(os.path.join(SRC, 'common_test.inc')).rstrip('\n').replace('{PU}', pre.upper()).replace('{LABEL}', name)
    cd = read(os.path.join(SRC, 'common_demo.inc')).rstrip('\n').replace('{P}', pre)
    rest = rest.replace('{COMMON_TEST}', ct).replace('{COMMON_DEMO}', cd).replace('{SETUP_CONSOLE}', SETUP_CONSOLE)
    assert '{' + 'COMMON' not in rest
    parts = [head, '', START] + [body(d) for d in deps] + [body(name), END, '', '', rest.lstrip('\n')]
    text = '\n'.join(parts)
    if not text.endswith('\n'):
        text += '\n'
    os.makedirs(OUT, exist_ok=True)
    out = os.path.join(OUT, num + '_' + name + '.cpp')
    with open(out, 'w', encoding='utf-8-sig', newline='\n') as f:
        f.write(text)
    return out


if __name__ == '__main__':
    names = sys.argv[1:] or [n for n in FILES if os.path.exists(os.path.join(SRC, n + '.body.hpp'))]
    for n in names:
        print(build(n))
