"""In E6, WHERE does simple spectrum fail?  Two questions:
 (1) is failure detected at the trivalent node alone (the analogue of
     thm:lemA, the reduction to t, which in type D is proved)?
 (2) among the 831 heavy-at-node weights, is the repeat always AT t?
usage: python3 ewhere.py [E6]
"""
import sys, collections
from rootsys import RootSystem, orbit_with_words
from wordtrace import trace
from socle import soc_sigma
T,N=(sys.argv[1][0],int(sys.argv[1][1:])) if len(sys.argv)>1 else ('E',6)
rs=RootSystem(T,N); t=[i for i in rs.I if len(rs.adj[i])>=3][0]
red=0; tab=collections.Counter(); heavy_rep=collections.Counter()
KR=[k for k in rs.I] if T!="D" else [k for k in rs.I if 2<=k<=N-2]
for k in KR:
    lam=rs.fundamental(k)
    for mu,w in orbit_with_words(rs,lam).items():
        if not w: continue
        tr=trace(rs,w,lam)
        _,rows,dv=soc_sigma(rs,k,w)
        # per-vertex multiplicity-freeness
        bad={i for d,row in rows.items() for i,c in row.items() if c>1}
        simple = not bad
        simple_at_t = t not in bad
        if simple != simple_at_t: red+=1
        tab[(simple, simple_at_t)]+=1
        if any(i==t and c>=2 for _,i,c,_ in tr):
            heavy_rep['heavy total']+=1
            heavy_rep['repeat occurs at t' if t in bad else
                      ('repeat elsewhere only' if bad else 'NO repeat (simple!)')]+=1
print(f"{T}{N}: trivalent node {t}")
print(f"(1) reduction to the node: 'simple' vs 'simple at {t}' disagree in {red} cases")
for kk in sorted(tab,key=str): print(f"      simple={kk[0]!s:5} simple_at_t={kk[1]!s:5} -> {tab[kk]}")
print("(2) among heavy-at-node weights:")
for kk in sorted(heavy_rep): print(f"      {kk:24} {heavy_rep[kk]}")
