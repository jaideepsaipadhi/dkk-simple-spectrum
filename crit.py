"""
REFINED CRITERION:
    NOT simple spectrum  <=>  c_max >= 2  AND  v_t >= 3   (t = trivalent node)
equivalently  simple  <=>  every c = 1  OR  v_t <= 2.
"""
from rootsys import RootSystem, orbit_with_words
from wordtrace import trace, v_from_trace
from dfs3 import solve

def run(typ, ns, cap=400):
    tot = exc = 0; bad = []
    for n in ns:
        rs = RootSystem(typ, n)
        tri = [i for i in rs.I if len(rs.adj[i]) >= 3]
        if not tri: continue
        t = tri[0]
        for k in rs.I:
            orb = orbit_with_words(rs, rs.fundamental(k))
            if len(orb) > cap: continue
            for mu, word in orb.items():
                if not word: continue
                tr = trace(rs, word, rs.fundamental(k))
                v = v_from_trace(rs, tr)
                cmax = max(c for _, _, c, _ in tr)
                r = solve(rs, k, v)
                if not r: continue
                rows = list(r.values())[0]
                ss = all(c <= 1 for x in rows.values() for c in x.values())
                pred = (cmax == 1) or (v[t-1] <= 2)
                tot += 1
                if ss != pred:
                    exc += 1
                    if len(bad) < 6:
                        bad.append((rs.name(), k, word, v, cmax, ss, pred))
    return tot, exc, bad

for typ, ns in [("D", range(4, 9)), ("E", [6, 7])]:
    tot, exc, bad = run(typ, ns)
    print(f"type {typ}: {tot} elements, {exc} exceptions to the refined criterion")
    for nm, k, w, v, cm, ss, pr in bad:
        print(f"   {nm} om{k} w={''.join('s'+str(i) for i in w)} v={v} "
              f"cmax={cm} simple={ss} pred={pr}")
