"""Validate every count asserted in appendixE.
usage: python3 validate_counts.py E6|E7|E8
"""
import sys, collections
from rootsys import RootSystem, orbit_with_words
from wordtrace import trace
from socle import soc_sigma
T_,N=sys.argv[1][0],int(sys.argv[1][1:])
rs=RootSystem(T_,N); t=[i for i in rs.I if len(rs.adj[i])>=3][0]
orbit={}; tot=0; has3=0; has4plus=0; max4=collections.Counter()
ones2_fail=0; nodeheavy=0; nodeheavy_simple=0; mx=0; mxsimple=0; nine=0
for k in rs.I:
    lam=rs.fundamental(k); O=orbit_with_words(rs,lam)
    orbit[k]=len(O)-1          # exclude the trivial (identity) weight
    for mu,w in O.items():
        if not w: continue
        tot+=1
        tr=trace(rs,w,lam); cs=[c for _,i,c,_ in tr]
        m=max(cs); mx=max(mx,m)
        _,rows,dv=soc_sigma(rs,k,w)
        simple=all(x<=1 for r in rows.values() for x in r.values())
        if simple: mxsimple=max(mxsimple,m)
        if m>=3: has3+=1
        if m>=4: has4plus+=1; max4[(k,m)]+=1
        if m==3 and simple: nine+=1
        if cs[:-1]==[1]*(len(cs)-1) and cs[-1]==2 and not simple: ones2_fail+=1
        if any(i==t and c>=2 for _,i,c,_ in tr):
            nodeheavy+=1; nodeheavy_simple+=simple
print(f"{T_}{N}:")
print(f"  orbit sizes (nontrivial) by kappa: {[orbit[k] for k in rs.I]}")
print(f"  total nontrivial extremal weights: {tot}")
print(f"  max c_p overall: {mx};  under simple spectrum: {mxsimple}")
print(f"  with a c-trace entry 3 (max c >= 3): {has3}")
print(f"  with max c >= 4: {has4plus}   by (kappa, maxc): {dict(max4)}")
print(f"  simple AND max c = 3 (the 'nine'): {nine}")
print(f"  c-trace (1,...,1,2) failing simple spectrum: {ones2_fail}")
print(f"  c-trace > 1 at the node: {nodeheavy}, of which simple: {nodeheavy_simple}")
