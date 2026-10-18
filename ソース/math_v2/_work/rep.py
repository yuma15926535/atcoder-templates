#!/usr/bin/env python3
"""rep.py <file> : stdin に「<<<OLD\n...\n===\n...\n>>>」の組を並べて置換（各 OLD はちょうど 1 回現れること）"""
import sys, re
p = sys.argv[1]
s = open(p, encoding='utf-8').read()
src = sys.stdin.read()
pairs = re.findall(r'<<<OLD\n(.*?)\n===\n(.*?)\n>>>', src, re.S)
assert pairs, 'no pairs'
for old, new in pairs:
    c = s.count(old)
    if c != 1:
        sys.exit(f'{p}: OLD が {c} 回: {old[:80]!r}')
    s = s.replace(old, new)
open(p, 'w', encoding='utf-8', newline='\n').write(s)
print(f'{p}: {len(pairs)} 箇所を置換')
