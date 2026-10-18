#!/usr/bin/env python3
"""Main_base（ACL 込み）に、グラフ・グリッド・幾何・数学・文字列の第二世代の本体を全部貼ってコンパイル・実行する（名前の衝突の確認。作業用）
使い方: check_all.py <v2 の親フォルダ> <文字列の out フォルダ> <ACL のフォルダ> <solv の中身のファイル>"""
import os, sys, glob, subprocess
START = 'ここからテンプレート本体（コピー範囲）'
END = 'ここまでテンプレート本体'
root, sout, acl, solvf = sys.argv[1:5]
def bodies(folder):
    res = []
    for f in sorted(glob.glob(os.path.join(folder, '*.cpp'))):
        t = open(f, encoding='utf-8-sig').read()
        a, b = t.find(START), t.find(END)
        if a < 0 or b < 0: continue
        a = t.index('\n', a) + 1
        b = t.rindex('\n', 0, b)
        res.append((os.path.basename(f), t[a:b]))
    return res
groups = [bodies(os.path.join(root, d)) for d in ('graph_grid_v2', 'geometry_v2', 'math_v2')] + [bodies(sout)]
print('本数:', [len(g) for g in groups])
main = open(os.path.join(root, 'Main_base.cpp'), encoding='utf-8-sig').read()
solv = open(solvf, encoding='utf-8').read()
k = main.index('void solv(){')
kend = main.index('int main(){')
ok = True
for order_name, order in (('順', groups), ('逆', [list(reversed(g)) for g in reversed(groups)])):
    code = main[:k] + '\n'.join(b for g in order for _, b in g) + '\n' + solv + '\n' + main[kend:]
    fn = f'/tmp/all_{order_name}.cpp'
    open(fn, 'w', encoding='utf-8').write(code)
    for cxx in ('g++', 'clang++'):
        for std in ('gnu++2b', 'c++17'):
            for local in ('', '-DLOCAL'):
                exe = '/tmp/all_exe'
                cmd = [cxx, f'-std={std}', '-O1', '-I', acl, '-o', exe, fn] + ([local] if local else [])
                r = subprocess.run(cmd, capture_output=True, text=True)
                if r.returncode:
                    ok = False
                    print('NG', order_name, cxx, std, local, r.stderr[:3000]); continue
                r2 = subprocess.run([exe], capture_output=True, text=True, timeout=120)
                print(order_name, cxx, std, local or 'nolocal', 'exit', r2.returncode, r2.stdout.strip().splitlines()[-1:] if r2.stdout.strip() else '')
                if r2.returncode: ok = False; print(r2.stderr[:1000])
print('ALL OK' if ok else 'FAILED')
