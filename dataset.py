"""Build the labelled dataset once so hypothesis search is fast."""
import json, time
from rootsys import RootSystem, orbit_with_words
from wordtrace import trace, v_from_trace
from dfs3 import solve

DATA = []
t0 = time.time()
for typ, ns in [("D", range(4, 9)), ("E", [6, 7]), ("A", range(3, 7))]:
    for n in ns:
        rs = RootSystem(typ, n)
        for k in rs.I:
            orb = orbit_with_words(rs, rs.fundamental(k))
            if len(orb) > 400: continue
            for mu, word in orb.items():
                if not word: continue
                tr = trace(rs, word, rs.fundamental(k))
                v = v_from_trace(rs, tr)
                cs = [c for _, _, c, _ in tr]
                r = solve(rs, k, v)
                if not r: continue
                rows = list(r.values())[0]
                ss = all(c <= 1 for x in rows.values() for c in x.values())
                DATA.append(dict(typ=typ, n=n, k=k, word=list(word),
                                 v=[int(x) for x in v], cmax=max(cs),
                                 ndef=sum(1 for c in cs if c >= 2),
                                 ell=len(word), h=max(rows), simple=ss))
json.dump(DATA, open("data.json", "w"))
print(f"{len(DATA)} labelled elements in {time.time()-t0:.0f}s")
from collections import Counter
print(Counter((d['typ'], d['simple']) for d in DATA))
