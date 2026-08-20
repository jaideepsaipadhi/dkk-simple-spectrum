"""
Local inequality at the trivalent node.

The module V is generated in degree 0, so for every d the map
    Z : (+)_{j~t} V_j[d+1]  ->  V_t[d+2]
is surjective, and the preprojective relation at t says Z o Y = 0 where
    Y : V_t[d]  ->  (+)_{j~t} V_j[d+1].
Hence  im(Y) <= ker(Z)  and

    dim V_t[d+2]  <=  SUM_{j~t} dim V_j[d+1]  -  dim V_t[d]  +  dim ker(Y).

If Y is injective this is the clean inequality

    dim V_t[d] + dim V_t[d+2]  <=  SUM_{j~t} dim V_j[d+1].          (*)

Test (*) on every element in the dataset.
"""
import json
from rootsys import RootSystem
from dfs3 import solve
ALL = json.load(open("data.json"))
RS={}
def rs_of(d):
    kk=(d['typ'],d['n'])
    if kk not in RS: RS[kk]=RootSystem(*kk)
    return RS[kk]

viol = tot = 0
tight = 0
ex=[]
for d in ALL:
    rs=rs_of(d); tri=[i for i in rs.I if len(rs.adj[i])>=3]
    if not tri: continue
    t=tri[0]
    rows=list(solve(rs,d['k'],tuple(d['v'])).values())[0]
    def dim(dd,i): return rows.get(dd,{}).get(i,0)
    top=max(rows)
    for dd in range(-2, top+1):
        lhs = dim(dd,t) + dim(dd+2,t)
        rhs = sum(dim(dd+1,j) for j in rs.adj[t])
        if lhs or rhs:
            tot += 1
            if lhs > rhs:
                viol += 1
                if len(ex)<4: ex.append((d['typ'],d['n'],d['k'],d['v'],dd,lhs,rhs))
            elif lhs == rhs: tight += 1
print(f"inequality (*) checked at {tot} (element, degree) pairs")
print(f"   violations: {viol}")
print(f"   equalities: {tight}")
for x in ex: print("   ",x)
