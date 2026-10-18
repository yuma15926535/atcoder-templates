import glob
START='ここからテンプレート本体（コピー範囲）'; END='// ======================= ここまでテンプレート本体'
ours=[]
for f in sorted(glob.glob('/home/claude/sv2/_work/out/*.cpp')):
    t=open(f,encoding='utf-8-sig').read()
    ours.append(t[t.index('\n',t.index(START))+1:t.index(END)])
s=open('cmp3.cpp',encoding='utf-8').read()
key='#define rep(i,n) for(ll i=0;i<(ll)(n);++i)\n'
a=s.index(key)+len(key); b=s.index('\nnamespace old_rh {')
s=s[:a]+'\n'.join(ours)+s[b:]
open('cmp3.cpp','w',encoding='utf-8').write(s)
