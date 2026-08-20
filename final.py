import sys, time
from rootsys import RootSystem, orbit_with_words
from wordtrace import trace, v_from_trace
from dfs3 import solve

def arms(rs, t):
    """arm parameters (p,q,r) at a trivalent node: arm of length L -> L+1"""
    out = []
    for nb in sorted(rs.adj[t]):
        seen = {t}; cur = nb; L = 1
        while True:
            nxt = [x for x in rs.adj[cur] if x not in seen]
            seen.add(cur)
            if len(nxt) != 1: break
            cur = nxt[0]; L += 1
        out.append(L + 1)
    return sorted(out)

CASES = [("D",4),("D",5),("D",6),("D",7),("D",8),("E",6),("E",7)]
for typ, n in CASES:
    rs = RootSystem(typ, n)
    tri = [i for i in rs.I if len(rs.adj[i]) >= 3]
    t0 = tri[0]
    ar = arms(rs, t0)
    print(f"\n=== {rs.name()}  trivalent node {t0}, arms {tuple(ar)}, "
          f"2nd smallest = {ar[1]}")
    for k in rs.I:
        orb = orbit_with_words(rs, rs.fundamental(k))
        if len(orb) > 400: continue
        S, NS = set(), set()
        unres = 0; t = time.time()
        for mu, word in orb.items():
            if not word: continue
            v = v_from_trace(rs, trace(rs, word, rs.fundamental(k)))
            r = solve(rs, k, v)
            if not r: unres += 1; continue
            rows = list(r.values())[0]
            ss = all(c <= 1 for x in rows.values() for c in x.values())
            (S if ss else NS).add(v[t0-1])
        ov = S & NS
        thr = max(S) if (S and not ov) else None
        note = ("all simple" if not NS else
                (f"threshold v_t <= {thr}" if thr is not None else f"OVERLAP {sorted(ov)}"))
        flag = ""
        if thr is not None:
            flag = "  <-- matches 2nd-smallest" if thr == ar[1] else f"  <-- 2nd-smallest is {ar[1]}"
        print(f"  om{k}: |orbit|={len(orb):4d} unres={unres:2d}  "
              f"simple v_t={sorted(S)} not={sorted(NS)}  {note}{flag}  [{time.time()-t:.1f}s]")
