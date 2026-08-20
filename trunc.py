"""
LOCALITY DEPTH.

If the classification depends only on v_t together with a BOUNDED prefix of
each arm's v-sequence (plus which arm holds k), then in type D the criterion
can be proved by checking finitely many local configurations, independent of
rank.  Find the minimal prefix depth with zero ambiguity.
"""
import json
from collections import defaultdict
from rootsys import RootSystem

ALL = json.load(open("data.json")) + [json.loads(l) for l in open("extra.jsonl")]
RS = {}
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

def keyfn(d, depth, use_k=True):
    rs=rs_of(d); tri=[i for i in rs.I if len(rs.adj[i])>=3]
    t=tri[0]
    parts=[]
    for A in arms(rs,t):
        seq=tuple(d['v'][x-1] for x in A[:depth])
        parts.append(((d['k'] in A) if use_k else False, seq, min(len(A), depth+1)))
    return (d['typ'], d['v'][t-1], tuple(sorted(parts)))

for typ in ("D","E"):
    sub=[d for d in ALL if d['typ']==typ and d['cmax']!=1
         and [i for i in rs_of(d).I if len(rs_of(d).adj[i])>=3]]
    print(f"\ntype {typ}: {len(sub)} elements with c_max>=2")
    for depth in range(1,6):
        m=defaultdict(set)
        for d in sub: m[keyfn(d,depth)].add(d['simple'])
        amb=[k for k,v in m.items() if len(v)>1]
        print(f"   prefix depth {depth}: {len(m):4d} keys, ambiguous {len(amb)}")
        if not amb:
            print(f"   --> DEPTH {depth} SUFFICES in type {typ}")
            break
