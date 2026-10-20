"""Type D's proof runs on (H): V is a top truncation at t and an UP-SET at each
leaf.  E6's node has one true leaf (vertex 2) and two长 arms.  Which vertices is
V an up-set at, and does that predict the column shape?

up-set at i: the occupied degrees at i are exactly the top ones, each slot FULL
   -- V_i[d] = I_kappa,i[d] for d >= D, and 0 for d < D.
usage: python3 eupset.py E6
"""
import sys, collections
from rootsys import RootSystem, orbit_with_words
from socle import soc_sigma, proj, nakayama
T_,N=(sys.argv[1][0],int(sys.argv[1][1:])) if len(sys.argv)>1 else ('E',6)
rs=RootSystem(T_,N); t=[i for i in rs.I if len(rs.adj[i])>=3][0]
def colprof(pr,i):
    return {d:pr.dim[(d,i)] for (d,j) in pr.dim if j==i and pr.dim[(d,i)]>0}
def gr(p,m):
    o={}
    for d in sorted(p,reverse=True):
        if m<=0: break
        x=min(p[d],m)
        if x: o[d]=x; m-=x
    return o,m
from upset import is_upset
tab=collections.Counter()
for k in rs.I:
    lam=rs.fundamental(k); pr=proj(rs,nakayama(rs)[k])
    PT=colprof(pr,t); ONES={d:1 for d in PT}
    for mu,w in orbit_with_words(rs,lam).items():
        if not w: continue
        _,rows,dv=soc_sigma(rs,k,w)
        col={d:r[t] for d,r in rows.items() if t in r}
        a,la=gr(PT,dv[t-1]); b,lb=gr(ONES,dv[t-1])
        shape=('full' if (col==a and la==0) else
               'strand' if (col==b and lb==0) else 'mixed')
        ups=tuple(i for i in rs.I
                  if is_upset({d:r[i] for d,r in rows.items() if i in r}, colprof(pr,i)))
        tab[(shape,ups)]+=1
print(f"{T_}{N} (node {t}): shape vs the set of vertices V is an up-set at")
for kk in sorted(tab,key=lambda x:(x[0],x[1])):
    print(f"   {kk[0]:7} upset_at={str(kk[1]):22} -> {tab[kk]}")
print("\nper-shape: which vertices are ALWAYS up-sets?")
for sh in ('full','strand','mixed'):
    sets=[kk[1] for kk in tab if kk[0]==sh]
    if not sets: continue
    always=set(rs.I)
    for s in sets: always &= set(s)
    print(f"   {sh:7}: always up-set at {sorted(always)}  ({sum(v for kk,v in tab.items() if kk[0]==sh)} cases)")
