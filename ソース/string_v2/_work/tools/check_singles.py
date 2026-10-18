#!/usr/bin/env python3
"""Main_base ＋ 単独テンプレート（binary_search・static_seq_toolkit・dynamic_unique_seq・segtree・lazy_segtree）＋ 文字列の本体を貼ってコンパイル（作業用）
使い方: check_singles.py <claude フォルダ> <ACL> <文字列の out の .cpp ...>"""
import subprocess, sys
root, acl, files = sys.argv[1], sys.argv[2], sys.argv[3:]
def rd(f): return open(f, encoding='utf-8-sig').read()
def lines(f, a, b): return '\n'.join(rd(f).split('\n')[a-1:b])
main = rd(f'{root}/Main_base.cpp'); k = main.index('void solv(){'); ke = main.index('int main(){')
bs = lines(f'{root}/binary_search_template/binary_search.cpp', 180, 1080)
t = rd(f'{root}/static_seq_template/static_seq_toolkit.cpp'); ss = t[t.index('#include <bits/stdc++.h>'):t.index('// ============================ StaticSeq ここまで')]
t = rd(f'{root}/dynamic_unique_seq_template/dynamic_unique_seq.cpp'); a = t.index('// ===================== テンプレート本体 ここから'); dus = t[a:t.index('// ===================== テンプレート本体 ここまで')]
seg = lines(f'{root}/segtree_templates/segtree.cpp', 20, 211); lseg = lines(f'{root}/segtree_templates/lazy_segtree.cpp', 21, 157)
ours = []
for f in files:
    o = rd(f); ours.append(o[o.index('\n', o.index('ここからテンプレート本体')) + 1:o.index('// ======================= ここまでテンプレート本体')])
solv = 'void solv(){ string S="abcba"; cout << is_palindrome(S) << "\\n"; }\n'
ok = True
for name, parts in (('順', [bs, ss, dus, seg, lseg] + ours), ('逆', ours[::-1] + [lseg, seg, dus, ss, bs])):
    open('/tmp/single.cpp', 'w').write(main[:k] + '\n'.join(parts) + '\n' + solv + main[ke:])
    for std in ('gnu++2b', 'c++17'):
        r = subprocess.run(['g++', f'-std={std}', '-O1', '-I', acl, '-o', '/tmp/single', '/tmp/single.cpp'], capture_output=True, text=True)
        print(name, std, 'compile', 'OK' if r.returncode == 0 else 'NG\n' + r.stderr[:2000])
        ok &= r.returncode == 0
print('ALL OK' if ok else 'FAILED')
