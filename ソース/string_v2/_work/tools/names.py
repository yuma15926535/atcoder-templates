import re,sys,glob,os
START='ここからテンプレート本体'
END='ここまでテンプレート本体'
def strip(s):
    s=re.sub(r'//[^\n]*','',s)
    s=re.sub(r'/\*.*?\*/','',s,flags=re.S)
    s=re.sub(r'"(\\.|[^"\\])*"','""',s)
    s=re.sub(r"'(\\.|[^'\\])*'","''",s)
    return s
def handle(st,out):
    st=re.sub(r'\s+',' ',st).strip()
    if not st: return
    for _ in range(3): st=re.sub(r'^template\s*<[^<>]*(<[^<>]*>[^<>]*)*>\s*','',st)
    st=re.sub(r'^\[\[[^\]]*\]\]\s*','',st)
    m=re.match(r'^(?:struct|class|union|enum(?: class)?)\s+(?:alignas\(\w+\)\s*)?(\w+)',st)
    if m: out.add(m.group(1)); return
    m=re.match(r'^using\s+(\w+)\s*=',st)
    if m: out.add(m.group(1)); return
    if st.startswith('using namespace') or st.startswith('static_assert') or st.startswith('extern') : return
    m=re.match(r'^(?:__extension__ )?typedef .*?(\w+)$',st)
    if m: out.add(m.group(1)); return
    if 'operator' in st.split('(')[0]: return
    m=re.search(r'(\w+)\s*\(',st)
    if m and '=' not in st.split('(')[0]:
        nm=m.group(1)
        if nm not in ('if','for','while','return','sizeof','decltype','alignas','noexcept'): out.add(nm); return
    m=re.match(r'^[\w:<>,\s\*&]*?\b(\w+)\s*(=|\[|\{|$)',st)
    if m: out.add(m.group(1))
def names_of(text):
    a=text.find(START); b=text.find(END)
    if a<0 or b<0: return set()
    body=text[text.find('\n',a)+1:b]
    macros=set(re.findall(r'^\s*#\s*define\s+(\w+)',body,re.M))
    body=re.sub(r'^\s*#[^\n]*','',strip(body),flags=re.M)
    out=set(macros); stack=[]; buf=''
    def at_top(): return all(x=='ns-pub' for x in stack)
    for c in body:
        if c=='{':
            last=buf.strip().split(';')[-1].strip()
            m=re.match(r'^(?:inline\s+)?namespace\s*(\w*)$',last)
            if m:
                nm=m.group(1); stack.append('ns-int' if (nm=='' or 'internal' in nm or nm=='std') else 'ns-pub')
            else:
                if at_top(): handle(last,out)
                stack.append('blk')
            buf=''; continue
        if c=='}':
            if stack: stack.pop()
            buf=''; continue
        if c==';':
            if at_top(): handle(buf.strip(),out)
            buf=''; continue
        if at_top(): buf+=c
    return out
if __name__=='__main__':
    res={}
    for folder in sys.argv[1:]:
        fs=[folder] if folder.endswith('.cpp') else sorted(glob.glob(os.path.join(folder,'*.cpp')))
        for f in fs:
            for nm in names_of(open(f,encoding='utf-8-sig').read()):
                res.setdefault(nm,set()).add(os.path.basename(f))
    for k in sorted(res): print(k, ' '.join(sorted(res[k]))[:150])
