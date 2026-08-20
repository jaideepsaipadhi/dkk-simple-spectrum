import json
from collections import defaultdict
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

tab=defaultdict(set)
for d in DATA:
    rs=rs_of(d); tri=[i for i in rs.I if len(rs.adj[i])>=3]
    if not tri or d['cmax']==1: continue
    t=tri[0]; ar=arms(rs,t)
    key=(d['v'][t-1], tuple(sorted((d['k'] in A, tuple(d['v'][x-1] for x in A)) for A in ar)))
    tab[key] |= {d['simple']}

print("v_t = 2 cases (arms shown as (contains_k, v-sequence)):")
for key in sorted(tab, key=lambda x:(x[0], x[1])):
    if key[0]!=2: continue
    lab=list(tab[key])[0]
    if not lab: print(f"   NOT simple: {key[1]}")
n2=[k for k in tab if k[0]==2]
print(f"   ({sum(1 for k in n2 if list(tab[k])[0])} simple / {len(n2)} keys with v_t=2)")

print("\nv_t = 3, arm lengths (1,2,2)  [E6-shape]:")
for key in sorted(tab, key=lambda x:str(x)):
    if key[0]!=3: continue
    L=sorted(len(a[1]) for a in key[1])
    if L!=[1,2,2]: continue
    print(f"   {'SIMPLE ' if list(tab[key])[0] else 'not    '} {key[1]}")
