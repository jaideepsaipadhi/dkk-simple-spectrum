"""Is simplicity determined by LOCAL data at the trivalent node?"""
import json
from collections import defaultdict
from rootsys import RootSystem
DATA = json.load(open("data.json"))
RS = {}
def rs_of(d):
    kk=(d['typ'],d['n'])
    if kk not in RS: RS[kk]=RootSystem(*kk)
    return RS[kk]

def arms(rs, t):
    out=[]
    for nb in sorted(rs.adj[t]):
        seen={t}; cur=nb; path=[nb]
        while True:
            nxt=[x for x in rs.adj[cur] if x not in seen]
            seen.add(cur)
            if len(nxt)!=1: break
            cur=nxt[0]; path.append(cur)
        out.append(path)
    return out

def test(name, keyfn):
    m = defaultdict(set)
    for d in DATA:
        rs = rs_of(d); tri=[i for i in rs.I if len(rs.adj[i])>=3]
        if not tri: continue
        if d['cmax'] == 1: continue          # these are always simple
        m[keyfn(d, rs, tri[0])].add(d['simple'])
    bad = {k: v for k, v in m.items() if len(v) > 1}
    print(f"{name:<52} keys {len(m):4d}, ambiguous {len(bad)}")
    for k in list(bad)[:4]: print(f"      ambiguous: {k}")
    return bad

test("v_t alone",
     lambda d, rs, t: (d['typ'], d['v'][t-1]))
test("v_t + sorted neighbour v",
     lambda d, rs, t: (d['typ'], d['v'][t-1], tuple(sorted(d['v'][j-1] for j in rs.adj[t]))))
test("v_t + neighbour v by arm length",
     lambda d, rs, t: (d['typ'], d['v'][t-1],
                       tuple(sorted((len(A), d['v'][A[0]-1]) for A in arms(rs, t)))))
test("full v-sequence along each arm",
     lambda d, rs, t: (d['v'][t-1],
                       tuple(sorted(tuple(d['v'][x-1] for x in A) for A in arms(rs, t)))))
