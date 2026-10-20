"""Cache (kappa, v, column shape, simple) for a type, so predicates can be
searched against E6 AND E7 jointly instead of fitted to one.
usage: python3 edump.py E6 e6shape.tsv
"""
import sys
from rootsys import RootSystem, orbit_with_words
from socle import soc_sigma, proj, nakayama
T_,N=sys.argv[1][0],int(sys.argv[1][1:]); out=sys.argv[2]
rs=RootSystem(T_,N); t=[i for i in rs.I if len(rs.adj[i])>=3][0]
def prof(k):
    pr=proj(rs,nakayama(rs)[k])
    return {d:pr.dim.get((d,t),0) for d in range(pr.top+1) if pr.dim.get((d,t),0)>0}
def gr(p,m):
    o={}
    for d in sorted(p,reverse=True):
        if m<=0: break
        x=min(p[d],m)
        if x: o[d]=x; m-=x
    return o,m
with open(out,'w') as f:
    f.write("kappa\tv\tshape\tsimple\tvt\n")
    for k in rs.I:
        lam=rs.fundamental(k); P=prof(k); ONES={d:1 for d in P}
        for mu,w in orbit_with_words(rs,lam).items():
            if not w: continue
            _,rows,dv=soc_sigma(rs,k,w)
            col={d:r[t] for d,r in rows.items() if t in r}
            a,la=gr(P,dv[t-1]); b,lb=gr(ONES,dv[t-1])
            shape=('full' if (col==a and la==0) else
                   'strand' if (col==b and lb==0) else 'mixed')
            s=all(c<=1 for r in rows.values() for c in r.values())
            f.write(f"{k}\t{'|'.join(map(str,dv))}\t{shape}\t{int(s)}\t{dv[t-1]}\n")
print(f"{T_}{N} -> {out}")
