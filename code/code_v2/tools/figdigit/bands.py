import sys, re, subprocess, numpy as np
import os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from pdfpaths import paths
def labels(fn):
    out=subprocess.run(['pdftotext','-bbox',fn,'-'],capture_output=True,text=True).stdout
    L=[]
    for m in re.finditer(r'<word xMin="([\d.]+)" yMin="([\d.]+)" xMax="([\d.]+)" yMax="([\d.]+)">(.*?)</word>',out):
        x0,y0,x1,y1=map(float,m.groups()[:4]); L.append(((x0+x1)/2+57,(y0+y1)/2+12,m.group(5)))  # -> stream coords
    return L
def bands(fn):
    P=paths(fn); labs=labels(fn)
    ticks=[p['pts'][0] for p in P if p['op']=='S' and p['scol']==(0.0,0.0,0.0) and len(p['pts'])==1 and len(p['pts'][0])==2]
    vt=[(a[0],a[1],b[1]) for a,b in ticks if abs(a[0]-b[0])<1e-3 and 7.5<abs(a[1]-b[1])<8.6]   # major vertical ticks (x axis)
    ht=[(a[1],a[0],b[0]) for a,b in ticks if abs(a[1]-b[1])<1e-3 and 7.5<abs(a[0]-b[0])<8.6]   # major horizontal ticks (y axis)
    edges=[p['pts'][0] for p in P if p['op']=='S' and p['scol']==(1.0,0.7,0.7)]
    res=[]
    for e in edges:
        e=np.array(e); xm=e[:,0].mean(); ym=e[:,1].mean()
        res.append((xm,ym,e))
    return vt,ht,res,labs
if __name__=='__main__':
    vt,ht,res,labs=bands(sys.argv[1])
    print('vertical major ticks', len(vt)); print(sorted(set((round(x,1),round(min(y0,y1),1)) for x,y0,y1 in vt))[:40])
    print('horizontal major ticks', len(ht)); print(sorted(set((round(min(x0,x1),1),round(y,1)) for y,x0,x1 in ht))[:60])
    for xm,ym,e in res: print('edge', round(xm,1), round(ym,1), len(e), e[:,0].min(), e[:,0].max())

def panels(fn, energies):
    """returns dict E -> (theta_grid, lo, hi, tickvals); energies row-major (top-left first)"""
    vt,ht,res,labs=bands(fn)
    xs=sorted(set(round(e[:,0].min(),2) for _,_,e in res)); xe=sorted(set(round(e[:,0].max(),2) for _,_,e in res))
    cols=list(zip(xs,xe))
    ys=sorted(set(round(min(y0,y1),1) for x,y0,y1 in vt))   # top/bottom frame tick bases per row
    ys=sorted(set(ys))
    rows=[(ys[2*i],ys[2*i+1]+9) for i in range(len(ys)//2)]
    nums=[(x,y,float(t.replace('\u2212','-'))) for x,y,t in labs if re.match(r'^-?[\d.]+$',t.replace('\u2212','-'))]
    out={}; k=0
    for (ylo,yhi) in rows:
        g=sorted(set(round(y,2) for y,x0,x1 in ht if ylo-1<y<yhi+1))
        vals=[]
        for y in g:
            cand=[(abs(ly-y),v) for lx,ly,v in nums if lx<120 and abs(ly-y)<8]
            vals.append(min(cand)[1] if cand else np.nan)
        g=np.array(g); vals=np.array(vals); ok=~np.isnan(vals)
        cf=np.polyfit(g[ok],vals[ok],1)
        resid=np.max(np.abs(np.polyval(cf,g[ok])-vals[ok]))
        for (x0,x1) in cols:
            E=energies[k]; k+=1
            ed=[e for xm,ym,e in res if abs(e[:,0].min()-x0)<0.5 and ylo-5<ym<yhi+5]
            th=np.linspace(0.5,179.5,180)
            curves=[]
            for e in ed:
                tx=(e[:,0]-x0)/(x1-x0)*180.0; vy=np.polyval(cf,e[:,1])
                o=np.argsort(tx); curves.append(np.interp(th,tx[o],vy[o]))
            curves=np.array(curves)
            out[E]=(th,curves.min(0),curves.max(0),(vals,resid,len(ed)))
    return out
