"""T(kappa) in closed form from the colour-t profile of I_kappa:
   a = number of leading multiplicity-1 slots from the top degree;
   T = a+1 if the next slot has multiplicity >= 2, else a.
Printed for every type with a branch node, E8 included -- E8's orbits are far
too large to enumerate, but its profiles are not.
usage: python3 ethresh.py
"""
from rootsys import RootSystem
from socle import proj, nakayama
def prof_of(rs,k,t):
    pr=proj(rs,nakayama(rs)[k])
    return [pr.dim[(d,t)] for d in sorted((d for (d,i) in pr.dim if i==t), reverse=True)
            if pr.dim[(d,t)]>0]
def greedy_T(p):
    best=0
    for m in range(sum(p)+1):
        rem=m; out=[]
        for c in p:
            if rem<=0: break
            out.append(min(c,rem)); rem-=min(c,rem)
        if rem==0 and all(x<=1 for x in out): best=m
        else: break
    return best
def closed_T(p):
    a=0
    while a<len(p) and p[a]==1: a+=1
    return a+1 if a<len(p) and p[a]>=2 else a
for T,N in (('D',5),('D',6),('D',7),('D',8),('E',6),('E',7),('E',8)):
    rs=RootSystem(T,N); t=[i for i in rs.I if len(rs.adj[i])>=3][0]
    print(f"{T}{N} (node {t}):")
    for k in rs.I:
        p=prof_of(rs,k,t); g,c=greedy_T(p),closed_T(p)
        print(f"   kappa={k}: profile {p}  T={g}" + ("" if g==c else f"  CLOSED FORM DISAGREES ({c})"))
