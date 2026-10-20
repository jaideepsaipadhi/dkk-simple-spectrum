"""Is the colour-t column of V always ONE of two truncations?

  (a) greedy top truncation of the full colour-t profile  (type D: prop:T), or
  (b) greedy top truncation of the ALL-ONES profile on the same degrees
      -- one element per degree, which is what a single tau-strand would give.

If so, 'spread vs fill' is a choice between two explicit shapes, and simple
spectrum is exactly 'the all-ones shape is the one that occurs, when v_t
exceeds what (a) allows'.
usage: python3 estrand.py E6 [E7]
"""
import sys, collections
from rootsys import RootSystem, orbit_with_words
from socle import soc_sigma, proj, nakayama
for arg in (sys.argv[1:] or ['E6']):
    T_,N=arg[0],int(arg[1:]); rs=RootSystem(T_,N)
    t=[i for i in rs.I if len(rs.adj[i])>=3][0]
    def prof(k):
        pr=proj(rs,nakayama(rs)[k])
        return {d:pr.dim.get((d,t),0) for d in range(pr.top+1) if pr.dim.get((d,t),0)>0}
    def gr(p,m):
        out={}
        for d in sorted(p,reverse=True):
            if m<=0: break
            x=min(p[d],m)
            if x: out[d]=x; m-=x
        return (out,m)
    tab=collections.Counter(); ex={}
    for k in rs.I:
        lam=rs.fundamental(k); P=prof(k); ONES={d:1 for d in P}
        for mu,w in orbit_with_words(rs,lam).items():
            if not w: continue
            _,rows,dv=soc_sigma(rs,k,w)
            col={d:r[t] for d,r in rows.items() if t in r}
            a,la=gr(P,dv[t-1]); b,lb=gr(ONES,dv[t-1])
            simple=all(c<=1 for r in rows.values() for c in r.values())
            lab=('full' if (col==a and la==0) else
                 'strand' if (col==b and lb==0) else 'NEITHER')
            tab[(lab,simple)]+=1
            ex.setdefault((lab,simple),(k,dv,dict(sorted(col.items(),reverse=True))))
    print(f"\n{T_}{N} (node {t}):")
    for kk in sorted(tab,key=str):
        print(f"   shape={kk[0]:8} simple={kk[1]!s:5} -> {tab[kk]:6}   e.g. kappa={ex[kk][0]} v={ex[kk][1]} col={ex[kk][2]}")
    n=sum(tab.values()); bad=sum(v for kk,v in tab.items() if kk[0]=='NEITHER')
    print(f"   every column is one of the two shapes: {n-bad}/{n}")
