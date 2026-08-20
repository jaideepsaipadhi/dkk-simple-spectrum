"""Extend the dataset with larger E orbits."""
import json, time
from rootsys import RootSystem, orbit_with_words
from wordtrace import trace, v_from_trace
from dfs3 import solve

DATA = json.load(open("data.json"))
have = {(d['typ'], d['n'], d['k']) for d in DATA}
NEW = []
t0 = time.time()
for typ, n, k in [("E",7,2), ("E",7,6), ("E",8,8), ("E",8,7)]:
    if (typ,n,k) in have: continue
    rs = RootSystem(typ, n); lam = rs.fundamental(k)
    orb = orbit_with_words(rs, lam)
    cnt = 0
    for mu, word in orb.items():
        if not word: continue
        tr = trace(rs, word, lam); v = v_from_trace(rs, tr)
        cs = [c for _,_,c,_ in tr]
        r = solve(rs, k, v)
        if not r: continue
        rows = list(r.values())[0]
        ss = all(c <= 1 for x in rows.values() for c in x.values())
        NEW.append(dict(typ=typ, n=n, k=k, word=list(word), v=[int(x) for x in v],
                        cmax=max(cs), ndef=sum(1 for c in cs if c>=2),
                        ell=len(word), h=max(rows), simple=ss))
        cnt += 1
    print(f"  {rs.name()} om{k}: {cnt} elements  ({time.time()-t0:.0f}s cumulative)")
json.dump(DATA + NEW, open("data2.json","w"))
print(f"total now {len(DATA)+len(NEW)}")
