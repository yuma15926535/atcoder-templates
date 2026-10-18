# 使い方: python3 run_file.py <problems ファイル> [問題名 ...]
# problems ファイルは problems.py と同じ形（P(name, files, solv, samples) を並べたもの）
import os, subprocess, sys, importlib.util
sys.path.insert(0, '/home/claude/work')
import paste_test as PT
path = sys.argv[1]; only = set(sys.argv[2:])
PROBLEMS = []
def P(name, files, solv, samples): PROBLEMS.append((name, files, solv, samples))
code = open(path, encoding='utf-8').read()
exec(compile(code.replace('PROBLEMS = []', ''), path, 'exec'), {'P': P, 'PROBLEMS': PROBLEMS})
base = open(PT.W + '/Main_base.cpp', encoding='utf-8-sig').read()
old = 'void solv(){\n    return;\n}'
out = '/tmp/prob_' + os.path.basename(path).replace('.py', '')
os.makedirs(out, exist_ok=True)
allok = True
for name, files, solv, samples in PROBLEMS:
    if only and name not in only: continue
    bodies = '\n'.join(PT.body(PT.W + '/math/out/' + f + '.cpp') for f in files)
    src = base.replace(old, bodies + '\n' + solv.strip('\n'))
    p = f'{out}/{name}.cpp'
    open(p, 'w', encoding='utf-8').write(src)
    r = subprocess.run(f'g++ -std=gnu++2b -O2 -I/home/claude/acl/ac-library -o {out}/{name} {p}', shell=True, capture_output=True, text=True)
    if r.returncode != 0:
        print(name, 'COMPILE FAIL'); print('\n'.join([l for l in r.stderr.split('\n') if 'error' in l][:8])); allok = False; continue
    res = []
    for inp, exp in samples:
        o = subprocess.run(f'{out}/{name}', shell=True, input=inp, capture_output=True, text=True, timeout=60)
        ok = o.stdout.split() == exp.split() and o.returncode == 0
        res.append('OK' if ok else 'NG')
        if not ok: allok = False; print('  got', repr(o.stdout[:300]), o.stderr[:300])
    print(f'{name:14s} {" ".join(res)}')
print('ALL OK' if allok else 'NG')
