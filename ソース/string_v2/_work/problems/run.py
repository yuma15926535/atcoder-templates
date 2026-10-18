#!/usr/bin/env python3
"""代表的な問題を「Main_base ＋ テンプレート本体だけ」で解いて、入出力例を確かめる（作業用）
使い方: python3 run.py <problems ファイル> <Main_base.cpp> <out フォルダ> <ACL フォルダ> [問題名 ...]
problems ファイルは P(name, files, solv, samples) を並べたもの（files は out の .cpp の名前。依存は同梱済みなので 1 つでよい）"""
import os, subprocess, sys
path, mainp, outdir, acl = sys.argv[1:5]
only = set(sys.argv[5:])
PROBLEMS = []
def P(name, files, solv, samples, check=None): PROBLEMS.append((name, files, solv, samples, check))  # check(入力, 出力, 期待) → bool（答えが複数あるとき）
exec(compile(open(path, encoding='utf-8').read(), path, 'exec'), {'P': P})
START, END = 'ここからテンプレート本体（コピー範囲）', '// ======================= ここまでテンプレート本体'
def body(f):
    t = open(f, encoding='utf-8-sig').read()
    return t[t.index('\n', t.index(START)) + 1:t.index(END)]
base = open(mainp, encoding='utf-8-sig').read()
old = 'void solv(){\n    return;\n}'
assert old in base
tmp = '/tmp/prob_' + os.path.basename(path).replace('.py', '')
os.makedirs(tmp, exist_ok=True)
allok = True
for name, files, solv, samples, check in PROBLEMS:
    if only and name not in only: continue
    src = base.replace(old, '\n'.join(body(os.path.join(outdir, f + '.cpp')) for f in files) + '\n' + solv.strip('\n'))
    p = f'{tmp}/{name}.cpp'
    open(p, 'w', encoding='utf-8').write(src)
    r = subprocess.run(['g++', '-std=gnu++2b', '-O2', '-I', acl, '-o', f'{tmp}/{name}', p], capture_output=True, text=True)
    if r.returncode != 0:
        print(name, 'COMPILE FAIL'); print('\n'.join([l for l in r.stderr.split('\n') if 'error' in l][:8])); allok = False; continue
    res = []
    for inp, exp in samples:
        o = subprocess.run([f'{tmp}/{name}'], input=inp, capture_output=True, text=True, timeout=60)
        ok = o.returncode == 0 and (check(inp, o.stdout, exp) if check else o.stdout.split() == exp.split())
        res.append('OK' if ok else 'NG')
        if not ok: allok = False; print('  got', repr(o.stdout[:300]), o.stderr[:300])
    print(f'{name:14s} {" ".join(res)}')
print('ALL OK' if allok else 'NG')
