import re, zlib, sys, numpy as np
def streams(fn):
    data=open(fn,'rb').read()
    out=[]
    for s in re.findall(rb'stream\r?\n(.*?)\r?\nendstream', data, re.S):
        try: out.append(zlib.decompress(s).decode('latin1'))
        except: pass
    return out
def paths(fn):
    """return list of dicts: op ('S' stroke or 'f' fill), color, lw, dash, pts (list of subpaths)"""
    res=[]
    for st in streams(fn):
        toks=st.split()
        stack=[]; cur=[]; sub=[]; scol=(0,0,0); fcol=(0,0,0); lw=1; dash=''
        gs=[]
        i=0
        while i<len(toks):
            tk=toks[i]
            if re.match(r'^-?[\d.]+$',tk): stack.append(float(tk)); i+=1; continue
            try:
              pass
            except Exception: pass
            if tk=='w' and not stack: i+=1; continue
            if tk in ('RG','rg') and len(stack)<3: stack=[]; i+=1; continue
            if tk in ('m','l','c') and len(stack)<2: stack=[]; i+=1; continue
            if tk=='m':
                if sub: cur.append(sub)
                sub=[(stack[-2],stack[-1])]
            elif tk=='l': sub.append((stack[-2],stack[-1]))
            elif tk=='c': sub.append((stack[-2],stack[-1]))
            elif tk=='re':
                x,y,w,h=stack[-4:]; cur.append([(x,y),(x+w,y),(x+w,y+h),(x,y+h)])
            elif tk=='h': pass
            elif tk in ('S','s','f','f*','B','b','n','W','W*'):
                if sub: cur.append(sub); sub=[]
                if tk in ('S','s','f','f*','B','b'):
                    res.append(dict(op=tk,scol=scol,fcol=fcol,lw=lw,dash=dash,pts=cur))
                if tk not in ('W','W*'): cur=[]
            elif tk=='RG': scol=tuple(stack[-3:])
            elif tk=='rg': fcol=tuple(stack[-3:])
            elif tk=='G': scol=(stack[-1],)*3
            elif tk=='g': fcol=(stack[-1],)*3
            elif tk=='w': lw=stack[-1]
            elif tk=='d':
                # dash array: tokens before like [ a b ] p d
                j=i-1; arr=[]
                while j>=0 and toks[j]!='[': j-=1
                dash=' '.join(toks[j:i])
            elif tk=='q': gs.append((scol,fcol,lw,dash))
            elif tk=='Q':
                if gs: scol,fcol,lw,dash=gs.pop()
            stack=[] if not re.match(r'^-?[\d.]+$',tk) else stack
            i+=1
    return res
if __name__=='__main__':
    P=paths(sys.argv[1])
    from collections import Counter
    c=Counter((p['op'],p['scol'] if p['op'] in('S','s') else p['fcol'],p['dash']) for p in P if sum(len(s) for s in p['pts'])>20)
    for k,v in c.most_common(40): print(k,v)
