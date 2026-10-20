"""Automated search for a simple function of v separating spread from fill
at kappa in {3,4,5}, v_t = 3 in E6.  158 weights, 14 spread.
usage: python3 esearch.py
"""
import collections, itertools
from rootsys import RootSystem, orbit_with_words
from socle import soc_sigma
rs=RootSystem('E',6); t=4
C=[]
for k in (3,4,5):
    lam=rs.fundamental(k)
    for mu,w in orbit_with_words(rs,lam).items():
        if not w: continue
        _,rows,dv=soc_sigma(rs,k,w)
        if dv[t-1]!=3: continue
        C.append((k,dv,all(c<=1 for r in rows.values() for c in r.values())))
print(f"{len(C)} weights, {sum(1 for c in C if c[2])} spread")
F={
 'v3+v5':            lambda k,v: v[2]+v[4],
 'v3+v5 mod 2':      lambda k,v: (v[2]+v[4])%2,
 '(v3 mod2, v5 mod2)':lambda k,v: (v[2]%2, v[4]%2),
 'v2':               lambda k,v: v[1],
 'v2 mod 2':         lambda k,v: v[1]%2,
 'v1+v6':            lambda k,v: v[0]+v[5],
 'max v':            lambda k,v: max(v),
 '(v3,v5) sorted':   lambda k,v: tuple(sorted((v[2],v[4]))),
 'v3+v5-v4':         lambda k,v: v[2]+v[4]-v[3],
 'all of v3,v5 even':lambda k,v: v[2]%2==0 and v[4]%2==0,
 '(kappa, v3+v5)':   lambda k,v: (k,v[2]+v[4]),
 '(kappa, maxv)':    lambda k,v: (k,max(v)),
 'sum v mod 2':      lambda k,v: sum(v)%2,
 'v1,v6 parity':     lambda k,v: (v[0]%2,v[5]%2),
 '(v2 mod2, v3+v5)': lambda k,v: (v[1]%2, v[2]+v[4]),
 'nbr sum v2+v3+v5': lambda k,v: v[1]+v[2]+v[4],
}
best=[]
for name,f in F.items():
    d=collections.defaultdict(lambda:[0,0])
    for k,v,s in C: d[f(k,v)][s]+=1
    mixed=sum(1 for vv in d.values() if vv[0] and vv[1])
    wrong=sum(min(vv) for vv in d.values())
    best.append((wrong,mixed,len(d),name,dict(d)))
best.sort()
for wrong,mixed,nd,name,d in best:
    flag='*** SEPARATES ***' if mixed==0 else f'{wrong} misclassified'
    print(f"\n{name}: {nd} values, {flag}")
    if mixed==0 or wrong<=14:
        for kk in sorted(d,key=str): print(f"    {kk!s:22} fill={d[kk][0]:3} spread={d[kk][1]:3}")
