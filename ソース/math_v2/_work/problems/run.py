import os, subprocess, sys
sys.path.insert(0, '/home/claude/work'); sys.path.insert(0, '/home/claude/work/problems')
import paste_test as PT
from problems import PROBLEMS
base = open(PT.W + '/Main_base.cpp', encoding='utf-8-sig').read()
old = 'void solv(){\n    return;\n}'
os.makedirs('/tmp/prob', exist_ok=True)
allok = True
for name, files, solv, samples in PROBLEMS:
    bodies = '\n'.join(PT.body(PT.W + '/math/out/' + f + '.cpp') for f in files)
    src = base.replace(old, bodies + '\n' + solv.strip('\n'))
    path = f'/tmp/prob/{name}.cpp'
    open(path, 'w', encoding='utf-8').write(src)
    r = subprocess.run(f'g++ -std=gnu++2b -O2 -I/home/claude/acl/ac-library -o /tmp/prob/{name} {path}', shell=True, capture_output=True, text=True)
    if r.returncode != 0:
        print(name, 'COMPILE FAIL'); print('\n'.join([l for l in r.stderr.split('\n') if 'error' in l][:5])); allok = False; continue
    res = []
    for i, (inp, exp) in enumerate(samples):
        o = subprocess.run(f'/tmp/prob/{name}', shell=True, input=inp, capture_output=True, text=True, timeout=60)
        ok = o.stdout.split() == exp.split() and o.returncode == 0
        res.append('OK' if ok else 'NG')
        if not ok: allok = False; print('  got', o.stdout[:200], o.stderr[:200])
    print(f'{name:12s} {" ".join(res)}')
print('ALL OK' if allok else 'NG')
