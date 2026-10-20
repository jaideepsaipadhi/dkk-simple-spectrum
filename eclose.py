"""The E6 endgame: at kappa in {3,4,5} with v_t = 3 = T+1, which weights have
simple spectrum?  Dump every candidate with the features that could separate
them, including the actual colour-t column against the greedy truncation.
usage: python3 eclose.py [E6]
"""
import sys, collections
from rootsys import RootSystem, orbit_with_words
from wordtrace import trace
from socle import soc_sigma, proj, nakayama
T,N=(sys.argv[1][0],int(sys.argv[1][1:])) if len(sys.argv)>1 else ('E',6)
rs=RootSystem(T,N); t=[i for i in rs.I if len(rs.adj[i])>=3][0]
def prof_of(k):
    pr=proj(rs,nakayama(rs)[k])
    return {d:pr.dim.get((d,t),0) for d in range(pr.top+1)}
def greedy(p,m):
    out={}
    for d in sorted(p,reverse=True):
        if m<=0: break
        take=min(p[d],m)
        if take: out[d]=take; m-=take
    return out
rows_out=[]
for k in rs.I:
    lam=rs.fundamental(k); P=prof_of(k)
    for mu,w in orbit_with_words(rs,lam).items():
        if not w: continue
        _,rows,dv=soc_sigma(rs,k,w)
        if dv[t-1]!=3 or k not in (3,4,5): continue
        simple=all(c<=1 for r in rows.values() for c in r.values())
        col={d:r[t] for d,r in rows.items() if t in r}
        g=greedy(P,3)
        rows_out.append(dict(ka=k,simple=simple,col=tuple(sorted(col.items(),reverse=True)),
                             greedy=tuple(sorted(g.items(),reverse=True)),
                             L=len(w),v=dv,mu=mu,
                             nb=(dv[2-1],dv[3-1],dv[5-1])))
print(f"{T}{N}: {len(rows_out)} weights with kappa in (3,4,5) and v_t=3; "
      f"{sum(r['simple'] for r in rows_out)} simple")
for key,name in [(lambda r:(r['ka'],r['col']),'(kappa, colour-t column)'),
                 (lambda r:r['col']==r['greedy'],'is V_t the greedy truncation?'),
                 (lambda r:(r['ka'],r['nb']),'(kappa, v at vertices 2,3,5)')]:
    d=collections.defaultdict(lambda:[0,0])
    for r in rows_out: d[key(r)][r['simple']]+=1
    mixed=[kk for kk,vv in d.items() if vv[0] and vv[1]]
    print(f"\n  by {name}: {len(d)} values, {len(mixed)} mixed"
          f" -> {'*** SEPARATES ***' if not mixed else 'mixed'}")
    for kk in sorted(d,key=str)[:16]:
        print(f"     {kk!s:46} not={d[kk][0]:3} simple={d[kk][1]:3}")
