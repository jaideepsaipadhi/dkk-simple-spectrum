"""Exhaustive verification of the type-D criterion."""
from rootsys import RootSystem, orbit_with_words
from wordtrace import trace, v_from_trace
from dfs3 import solve
tot = exc = 0
for n in range(4, 9):
    rs = RootSystem("D", n); t = n - 2
    for k in rs.I:
        orb = orbit_with_words(rs, rs.fundamental(k))
        if len(orb) > 400: continue
        for mu, word in orb.items():
            if not word: continue
            v = v_from_trace(rs, trace(rs, word, rs.fundamental(k)))
            r = solve(rs, k, v)
            if not r: continue
            rows = list(r.values())[0]
            ss = all(c <= 1 for x in rows.values() for c in x.values())
            tot += 1
            if ss != (v[t-1] <= 2):
                exc += 1
                if exc <= 5:
                    print(f"  EXCEPTION D{n} om{k} w={''.join('s'+str(i) for i in word)} "
                          f"v={v} simple={ss}")
print(f"\nType D, ranks 4-8, all fundamental weights (orbits <=400):")
print(f"  {tot} elements tested, {exc} exceptions to  'simple <=> v_t <= 2'")
