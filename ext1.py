import json, sys, time
from rootsys import RootSystem, orbit_with_words
from wordtrace import trace, v_from_trace
from dfs3 import solve
typ, n, k = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
rs = RootSystem(typ, n); lam = rs.fundamental(k)
out = open("extra.jsonl", "a")
t0=time.time(); cnt=0
for mu, word in orbit_with_words(rs, lam).items():
    if not word: continue
    tr = trace(rs, word, lam); v = v_from_trace(rs, tr)
    cs = [c for _,_,c,_ in tr]
    r = solve(rs, k, v)
    if not r: continue
    rows = list(r.values())[0]
    ss = all(c <= 1 for x in rows.values() for c in x.values())
    out.write(json.dumps(dict(typ=typ,n=n,k=k,word=list(word),v=[int(x) for x in v],
              cmax=max(cs), ell=len(word), h=max(rows), simple=ss))+"\n")
    cnt+=1
out.close()
print(f"{rs.name()} om{k}: {cnt} elements in {time.time()-t0:.0f}s")
