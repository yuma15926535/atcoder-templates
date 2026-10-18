#!/usr/bin/env python3
"""数学第二世代テンプレートの .cpp を組み立てる（作業用）。
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
    'int_common': ('00', [], 'intc'),
    'math_base': ('00', ['int_common'], 'mathb'),
    'modint': ('00', ['int_common', 'math_base'], 'mint'),
    'digits': ('01', ['int_common', 'math_base'], 'dig'),
    'factorize': ('02', ['int_common', 'math_base'], 'fac'),
    'sieve': ('02', ['int_common', 'math_base', 'factorize'], 'sieve'),
    'mod_equation': ('03', ['int_common', 'math_base'], 'meq'),
    'binomial': ('04', ['int_common', 'math_base', 'modint'], 'binom'),
    'matrix': ('05', ['int_common', 'math_base', 'modint'], 'mat'),
    'quotient': ('06', ['int_common', 'math_base'], 'quot'),
    'series': ('06', ['int_common', 'math_base', 'modint'], 'ser'),
    'fraction': ('08', ['int_common', 'math_base'], 'frac'),
    'bit': ('09', ['int_common', 'math_base'], 'bit'),
    # Phase 4（☆）
    'convolution': ('07', ['int_common', 'math_base', 'modint'], 'conv'),
    'fps': ('07', ['int_common', 'math_base', 'modint', 'convolution'], 'fps'),
    'subset_transform': ('07', ['int_common', 'math_base', 'modint'], 'subt'),
    'divisor_transform': ('02', ['int_common', 'math_base', 'factorize', 'sieve'], 'divt'),
    'prime_count': ('02', ['int_common', 'math_base'], 'pcnt'),
    'discrete_log': ('03', ['int_common', 'math_base', 'factorize'], 'dlog'),
    'combinatorics': ('04', ['int_common', 'math_base', 'modint', 'binomial', 'convolution'], 'comb'),
    'permutation': ('04', ['int_common', 'math_base', 'modint', 'binomial'], 'perm'),
    'linear_algebra': ('05', ['int_common', 'math_base', 'modint', 'matrix'], 'linalg'),
    'xor_basis': ('05', ['int_common', 'math_base'], 'xorb'),
    'linear_recurrence': ('06', ['int_common', 'math_base', 'modint', 'convolution'], 'linrec'),
    'bigint': ('08', ['int_common', 'math_base'], 'big'),
    'game': ('10', ['int_common', 'math_base'], 'game'),
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
