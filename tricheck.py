"""Is "c_p = 1 at every occurrence of the trivalent node" necessary for simple
spectrum?  Found in E6; tested here across type D, where the criterion r+2m<=2
is proved, so the verdict needs no recomputation of the theorem."""
import sys
from rootsys import RootSystem, orbit_with_words
from wordtrace import trace
from socle import soc_sigma
for typ in (('D',4),('D',5),('D',6),('D',7)):
    rs=RootSystem(*typ); t=[i for i in rs.I if len(rs.adj[i])>=3][0]
    tot=viol=simple_with_big=0
    for k in rs.I:
        lam=rs.fundamental(k)
        for mu,w in orbit_with_words(rs,lam).items():
            if not w: continue
            tr=trace(rs,w,lam)
            big_at_t=any(i==t and c>=2 for _,i,c,_ in tr)
            _,rows,_=soc_sigma(rs,k,w)
            simple=all(c<=1 for x in rows.values() for c in x.values())
            tot+=1
            if big_at_t: viol+=1; simple_with_big+=simple
    print(f"{typ[0]}{typ[1]}: trivalent node {t}; {tot} extremal weights, "
          f"{viol} with c_p>=2 at the node, of which simple: {simple_with_big}")
