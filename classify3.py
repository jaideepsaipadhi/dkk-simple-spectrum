import sys, time
from rootsys import RootSystem, orbit_with_words
from wordtrace import trace, v_from_trace
from dfs2 import solve

typ = sys.argv[1] if len(sys.argv) > 3 else "D"
n = int(sys.argv[2]) if len(sys.argv) > 3 else 5
k = int(sys.argv[3]) if len(sys.argv) > 3 else 2
dmax = int(sys.argv[4]) if len(sys.argv) > 4 else 14
rs = RootSystem(typ, n); lam = rs.fundamental(k)
tri = [i for i in rs.I if len(rs.adj[i]) >= 3]
orb = orbit_with_words(rs, lam)
simple_vt, notsimple_vt = set(), set()
res_n = unres = 0
t0 = time.time()
for mu, word in sorted(orb.items(), key=lambda t: len(t[1])):
    if not word: continue
    v = v_from_trace(rs, trace(rs, word, lam))
    r = solve(rs, k, v, dmax=dmax)
    if not r: unres += 1; continue
    res_n += 1
    rows = list(r.values())[0]
    ss = all(c <= 1 for x in rows.values() for c in x.values())
    for t in tri:
        (simple_vt if ss else notsimple_vt).add(v[t-1])
print(f"{rs.name()} om{k}: resolved {res_n}, unresolved {unres}, {time.time()-t0:.0f}s")
print(f"  trivalent node {tri}:  simple v_t = {sorted(simple_vt)}")
print(f"                        NOT simple v_t = {sorted(notsimple_vt)}")
ov = simple_vt & notsimple_vt
print(f"  OVERLAP: {sorted(ov) if ov else 'NONE -- clean threshold'}")
if not ov and simple_vt and notsimple_vt:
    print(f"  threshold: simple iff v_t <= {max(simple_vt)}")
