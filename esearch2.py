"""The last ambiguous class: (v3,v4,v5) = (2,3,2) in E6, kappa in {3,4,5}.
26 weights, 10 spread.  Search features on the arm ends v1, v2, v6.
usage: python3 esearch2.py
"""
import collections
from rootsys import RootSystem, orbit_with_words
from socle import soc_sigma
rs=RootSystem('E',6); t=4
C=[]
for k in (3,4,5):
    lam=rs.fundamental(k)
    for mu,w in orbit_with_words(rs,lam).items():
        if not w: continue
        _,rows,dv=soc_sigma(rs,k,w)
        if dv[t-1]!=3 or dv[2]!=2 or dv[4]!=2: continue
        C.append((k,dv,all(c<=1 for r in rows.values() for c in r.values())))
print(f"{len(C)} weights with (v3,v4,v5)=(2,3,2); {sum(1 for c in C if c[2])} spread\n")
for k,v,s in sorted(C,key=lambda r:(r[0],r[1])):
    print(f"  k={k} v={v} {'SPREAD' if s else 'fill  '}")
F={
 '(v1 mod2, v6 mod2)': lambda k,v:(v[0]%2,v[5]%2),
 '(kappa, v1 mod2, v6 mod2)': lambda k,v:(k,v[0]%2,v[5]%2),
 'v1+v6 mod 2':       lambda k,v:(v[0]+v[5])%2,
 'v2 mod 2':          lambda k,v:v[1]%2,
 '(v1+v6) ':          lambda k,v:v[0]+v[5],
 '(kappa, v2)':       lambda k,v:(k,v[1]),
 'v1*v6 parity pair sorted': lambda k,v:tuple(sorted((v[0]%2,v[5]%2))),
 '(v1 mod2,v2 mod2,v6 mod2)': lambda k,v:(v[0]%2,v[1]%2,v[5]%2),
}
print()
res=[]
for name,f in F.items():
    d=collections.defaultdict(lambda:[0,0])
    for k,v,s in C: d[f(k,v)][s]+=1
    wrong=sum(min(x) for x in d.values())
    res.append((wrong,name,dict(d)))
res.sort()
for wrong,name,d in res:
    print(f"{name}: {'*** SEPARATES ***' if wrong==0 else f'{wrong} misclassified'}")
    if wrong<=4:
        for kk in sorted(d,key=str): print(f"    {kk!s:24} fill={d[kk][0]:3} spread={d[kk][1]:3}")
