"""The exceptions to 'v_t <= T(kappa)': simple spectrum ABOVE the threshold.
These are exactly the weights where V_t is not the greedy top truncation, so
they are the whole remaining gap.
usage: python3 eexc.py [E6]
"""
import sys, collections
from rootsys import RootSystem, orbit_with_words
from wordtrace import trace
from socle import soc_sigma, proj, nakayama
T,N=(sys.argv[1][0],int(sys.argv[1][1:])) if len(sys.argv)>1 else ('E',6)
rs=RootSystem(T,N); t=[i for i in rs.I if len(rs.adj[i])>=3][0]
def prof_of(k):
    pr=proj(rs,nakayama(rs)[k])
    return [pr.dim[(d,t)] for d in sorted((d for (d,i) in pr.dim if i==t),reverse=True) if pr.dim[(d,t)]>0]
def Tof(p):
    a=0
    while a<len(p) and p[a]==1: a+=1
    return a+1 if a<len(p) and p[a]>=2 else a
th={k:Tof(prof_of(k)) for k in rs.I}
E=[]
for k in rs.I:
    lam=rs.fundamental(k)
    for mu,w in orbit_with_words(rs,lam).items():
        if not w: continue
        _,rows,dv=soc_sigma(rs,k,w)
        simple=all(c<=1 for r in rows.values() for c in r.values())
        if simple and dv[t-1]>th[k]:
            tr=trace(rs,w,lam)
            E.append(dict(ka=k,vt=dv[t-1],T=th[k],L=len(w),
                          nodec=tuple(c for _,i,c,_ in tr if i==t),v=dv))
print(f"{T}{N}: {len(E)} exceptions (simple above threshold)")
for key,name in [(lambda r:(r['ka'],r['vt'],r['T']),'(kappa, v_t, T)'),
                 (lambda r:r['nodec'],'c-values at the node'),
                 (lambda r:r['vt']-r['T'],'v_t - T')]:
    c=collections.Counter(key(r) for r in E)
    print(f"\n  by {name}:")
    for kk in sorted(c,key=str): print(f"     {kk!s:24} {c[kk]}")
