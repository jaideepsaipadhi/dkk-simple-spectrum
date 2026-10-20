"""The E6 failures of (*o): extremal V with c_t >= 1 that are NOT up-sets at a
neighbour j of t.  In type D prop:sharp says every failure has c_t = 0, so
these have no type-D analogue.  What are they?
usage: python3 estar2.py
"""
import collections
from rootsys import RootSystem, orbit_with_words
from socle import soc_sigma, proj, nakayama
from upset import is_upset
rs=RootSystem('E',6); t=4
def cp(pr,i): return {d:pr.dim[(d,i)] for (d,j) in pr.dim if j==i and pr.dim[(d,i)]>0}
def gr(p,m):
    o={}
    for d in sorted(p,reverse=True):
        if m<=0: break
        x=min(p[d],m)
        if x: o[d]=x; m-=x
    return o,m
bad=collections.defaultdict(list); shape_of={}
allct=collections.Counter()
for k in rs.I:
    lam=rs.fundamental(k); pr=proj(rs,nakayama(rs)[k])
    PT=cp(pr,t); ONES={d:1 for d in PT}
    for mu,w in orbit_with_words(rs,lam).items():
        if not w: continue
        _,rows,dv=soc_sigma(rs,k,w)
        ct=mu[t-1]
        col={d:r[t] for d,r in rows.items() if t in r}
        a,la=gr(PT,dv[t-1]); b,lb=gr(ONES,dv[t-1])
        sh=('full' if (col==a and la==0) else 'strand' if (col==b and lb==0) else 'mixed')
        for j in (2,3,5):
            u=is_upset({d:r[j] for d,r in rows.items() if j in r}, cp(pr,j))
            if not u:
                bad[j].append((k,ct,sh,dv))
                allct[(j,ct>=1,sh)]+=1
for j in (2,3,5):
    B=bad[j]; pos=[r for r in B if r[1]>=1]
    print(f"\nj={j}: {len(B)} non-up-sets total, {len(pos)} of them with c_t>=1")
    print(f"   by (c_t, shape) among ALL failures: "
          f"{dict(collections.Counter((r[1],r[2]) for r in B))}")
    if pos:
        print(f"   the c_t>=1 ones, by (kappa, c_t, shape): "
              f"{dict(collections.Counter((r[0],r[1],r[2]) for r in pos))}")
print("\n--- is 'non-full shape' the same set as 'fails (*o) somewhere with c_t>=1'? ---")
S1=set()
for j in (2,3,5):
    for k,ct,sh,dv in bad[j]:
        if ct>=1: S1.add((k,tuple(dv)))
print(f"   fails (*o) with c_t>=1: {len(S1)} weights; of these non-full: "
      f"{sum(1 for x in S1 if True)} (shapes: "
      f"{dict(collections.Counter(sh for j in (2,3,5) for k,ct,sh,dv in bad[j] if ct>=1))})")
