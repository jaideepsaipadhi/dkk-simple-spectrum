"""Cache the (*o) data: for every extremal V, the weight coordinates c_i =
<alpha_i^vee, mu>, the column shape at t, and whether V is an up-set at each
vertex.  Lets side conditions be searched against E6 AND E7 at once.
usage: python3 esdump.py E6 e6star.tsv
"""
import sys
from rootsys import RootSystem, orbit_with_words
from socle import soc_sigma, proj, nakayama
from upset import is_upset
T_,N=sys.argv[1][0],int(sys.argv[1][1:]); out=sys.argv[2]
rs=RootSystem(T_,N); t=[i for i in rs.I if len(rs.adj[i])>=3][0]
def cp(pr,i): return {d:pr.dim[(d,i)] for (d,j) in pr.dim if j==i and pr.dim[(d,i)]>0}
def gr(p,m):
    o={}
    for d in sorted(p,reverse=True):
        if m<=0: break
        x=min(p[d],m)
        if x: o[d]=x; m-=x
    return o,m
with open(out,'w') as f:
    f.write("kappa\tc\tv\tshape\tupsets\tsimple\n")
    for k in rs.I:
        lam=rs.fundamental(k); pr=proj(rs,nakayama(rs)[k])
        PT=cp(pr,t); ONES={d:1 for d in PT}
        profs={i:cp(pr,i) for i in rs.I}
        for mu,w in orbit_with_words(rs,lam).items():
            if not w: continue
            _,rows,dv=soc_sigma(rs,k,w)
            col={d:r[t] for d,r in rows.items() if t in r}
            a,la=gr(PT,dv[t-1]); b,lb=gr(ONES,dv[t-1])
            sh=('full' if (col==a and la==0) else 'strand' if (col==b and lb==0) else 'mixed')
            ups=[i for i in rs.I
                 if is_upset({d:r[i] for d,r in rows.items() if i in r}, profs[i])]
            s=all(c<=1 for r in rows.values() for c in r.values())
            f.write(f"{k}\t{'|'.join(map(str,mu))}\t{'|'.join(map(str,dv))}\t{sh}"
                    f"\t{'|'.join(map(str,ups))}\t{int(s)}\n")
print(f"{T_}{N} -> {out}")
