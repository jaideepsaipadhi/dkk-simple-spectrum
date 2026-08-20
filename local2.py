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
def kinfo(rs,t,k,ar):
    if k==t: return ("t",0)
    for a,A in enumerate(ar):
        if k in A: return (len(A), A.index(k))
    return ("?",0)

def test(name, keyfn):
    m=defaultdict(set)
    for d in DATA:
        rs=rs_of(d); tri=[i for i in rs.I if len(rs.adj[i])>=3]
        if not tri or d['cmax']==1: continue
        m[keyfn(d,rs,tri[0])].add(d['simple'])
    bad={k:v for k,v in m.items() if len(v)>1}
    print(f"{name:<58} keys {len(m):4d}  ambiguous {len(bad)}")
    for k in list(bad)[:5]: print(f"       {k}")
    return bad

test("v_t + arm v-sequences + k position",
     lambda d,rs,t: (d['v'][t-1],
                     tuple(sorted(tuple(d['v'][x-1] for x in A) for A in arms(rs,t))),
                     kinfo(rs,t,d['k'],arms(rs,t))))
test("v_t + arm v-seq (arms tagged by whether they contain k)",
     lambda d,rs,t: (d['v'][t-1],
                     tuple(sorted((d['k'] in A, tuple(d['v'][x-1] for x in A))
                                  for A in arms(rs,t)))))
