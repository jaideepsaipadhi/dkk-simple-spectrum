"""Where does a c-trace value >= 3 sit, in the simple-spectrum elements that
have one?  cor:nodeone (proved in type D) forbids c_p >= 2 at the trivalent
node under simple spectrum, so every such value should sit AWAY from it.
usage: python3 enine2.py E6 | E7
"""
import sys, collections
from rootsys import RootSystem, orbit_with_words
from wordtrace import trace
from socle import soc_sigma
T_,N=sys.argv[1][0],int(sys.argv[1][1:])
rs=RootSystem(T_,N); t=[i for i in rs.I if len(rs.adj[i])>=3][0]
hits=[]; maxsimple=0; maxall=0
for k in rs.I:
    lam=rs.fundamental(k)
    for mu,w in orbit_with_words(rs,lam).items():
        if not w: continue
        tr=trace(rs,w,lam); cs=[(i,c) for _,i,c,_ in tr]
        mx=max(c for _,c in cs); maxall=max(maxall,mx)
        _,rows,dv=soc_sigma(rs,k,w)
        if not all(x<=1 for r in rows.values() for x in r.values()): continue
        maxsimple=max(maxsimple,mx)
        if mx>=3:
            at=sorted({i for i,c in cs if c>=3}); node=[c for i,c in cs if i==t]
            hits.append((k,len(w),mx,at,node,dv))
print(f"{T_}{N} (node {t}): max c_p over ALL extremal weights = {maxall}; "
      f"over SIMPLE ones = {maxsimple}")
print(f"simple-spectrum weights with some c_p >= 3: {len(hits)}")
for k,L,mx,at,node,dv in hits:
    print(f"   kappa={k} len={L} max c={mx} at vertices {at}; c at node {t}: {node}; v={dv}")
bad=[h for h in hits if t in h[3]]
print(f"\nany with c_p >= 3 AT the node {t}: {len(bad)}")
