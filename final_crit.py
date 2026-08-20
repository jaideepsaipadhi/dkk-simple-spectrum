"""
CANDIDATE CRITERION (all ADE):

  simple spectrum  <=>  c_max = 1
                        OR  v_t <= 2
                        OR  some arm A at the trivalent node t with k not in A
                            carries the exact descending staircase
                            (v_t - 1, v_t - 2, ..., 1)

The "k not in A" clause is what makes type D fail at v_t = 3: there the only
arm long enough is the one containing the framing vertex.
"""
import json
from rootsys import RootSystem
DATA=json.load(open("data.json")); RS={}
def rs_of(d):
    kk=(d['typ'],d['n'])
    if kk not in RS: RS[kk]=RootSystem(*kk)
    return RS[kk]
def arms(rs,t):
    out=[]
    for nb in sorted(rs.adj[t]):
        seen={t};cur=nb;path=[nb]
        while True:
            nxt=[x for x in rs.adj[cur] if x not in seen]
            seen.add(cur)
            if len(nxt)!=1: break
            cur=nxt[0];path.append(cur)
        out.append(path)
    return out

def pred(d):
    rs=rs_of(d); tri=[i for i in rs.I if len(rs.adj[i])>=3]
    if d['cmax']==1 or not tri: return True
    t=tri[0]; v=d['v']; vt=v[t-1]
    if vt<=2: return True
    for A in arms(rs,t):
        if d['k'] in A: continue
        if len(A) < vt-1: continue
        if all(v[A[j]-1] == vt-1-j for j in range(vt-1)): return True
    return False

exc=0; bytyp={}
first=[]
for d in DATA:
    p=pred(d)
    if p!=d['simple']:
        exc+=1; bytyp[d['typ']]=bytyp.get(d['typ'],0)+1
        if len(first)<6:
            first.append((d['typ'],d['n'],d['k'],d['v'],d['cmax'],d['simple'],p))
print(f"CRITERION: {exc} exceptions out of {len(DATA)}   by type {bytyp}")
for x in first: print("   ",x)
