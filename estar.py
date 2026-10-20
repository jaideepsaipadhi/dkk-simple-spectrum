"""The (*o) analogue with three arms.

Type D (*o): an extremal M with c_t = <alpha_t^vee, lambda - dim M> >= 1 is an
up-set at the chain neighbour t-1.  Test, for every neighbour j of t, whether
   c_t >= 1  ==>  V is an up-set at j
holds; type D is the control, where it must hold at j = t-1.
usage: python3 estar.py D6 E6
"""
import sys, collections
from rootsys import RootSystem, orbit_with_words
from socle import soc_sigma, proj, nakayama
def colprof(pr,i): return {d:pr.dim[(d,i)] for (d,j) in pr.dim if j==i and pr.dim[(d,i)]>0}
from upset import is_upset
for arg in (sys.argv[1:] or ['D6','E6']):
    T_,N=arg[0],int(arg[1:]); rs=RootSystem(T_,N)
    t=[i for i in rs.I if len(rs.adj[i])>=3][0]
    nbrs=sorted(rs.adj[t])
    tab={j:collections.Counter() for j in nbrs}
    for k in rs.I:
        lam=rs.fundamental(k); pr=proj(rs,nakayama(rs)[k])
        for mu,w in orbit_with_words(rs,lam).items():
            if not w: continue
            _,rows,dv=soc_sigma(rs,k,w)
            ct=mu[t-1]              # <alpha_t^vee, mu> = c_t
            for j in nbrs:
                u=is_upset({d:r[j] for d,r in rows.items() if j in r}, colprof(pr,j))
                tab[j][(ct>=1,u)]+=1
    print(f"\n{T_}{N} (node {t}, neighbours {nbrs}):")
    for j in nbrs:
        c=tab[j]
        viol=c[(True,False)]
        tot=c[(True,False)]+c[(True,True)]
        print(f"   j={j}: among the {tot} extremal V with c_t>=1, "
              f"up-set at j in {c[(True,True)]}, FAILS in {viol}"
              + ("   <-- (*o) holds" if viol==0 else ""))
