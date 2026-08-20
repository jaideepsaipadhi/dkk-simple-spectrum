"""Where does the collision happen? Which vertex carries the repeated weight?"""
import json
from collections import Counter
from rootsys import RootSystem
from wordtrace import trace, v_from_trace
from dfs3 import solve

DATA = json.load(open("data.json"))
RS = {}
def rs_of(d):
    key=(d['typ'],d['n'])
    if key not in RS: RS[key]=RootSystem(*key)
    return RS[key]

cnt = Counter()
for d in DATA:
    if d['simple']: continue
    rs = rs_of(d); tri=[i for i in rs.I if len(rs.adj[i])>=3]
    t = tri[0] if tri else None
    rows = list(solve(rs, d['k'], tuple(d['v'])).values())[0]
    bad = {i for r in rows.values() for i, c in r.items() if c >= 2}
    cnt[(d['typ'], tuple(sorted('t' if i==t else ('nbr' if t and i in rs.adj[t] else 'far')
                                for i in bad)))] += 1
print("Which vertices carry a repeated weight, relative to the trivalent node:")
for kk, vv in cnt.most_common(12):
    print(f"   {kk}: {vv}")

print("\nMax multiplicity seen:")
mx = Counter()
for d in DATA:
    if d['simple']: continue
    rows = list(solve(rs_of(d), d['k'], tuple(d['v'])).values())[0]
    mx[max(c for r in rows.values() for c in r.values())] += 1
print("  ", dict(mx))
