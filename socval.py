"""Validate soc_sigma against dfs3.solve over whole orbits, including vertices NOT fixed
by the Nakayama permutation.  usage: python3 socval.py TYPE k [maxlen]"""
import sys, time
from rootsys import RootSystem, orbit_with_words
from wordtrace import trace, v_from_trace
from socle import soc_sigma, nakayama
from dfs3 import solve
typ = sys.argv[1]; T, N = typ[0], int(typ[1:]); k = int(sys.argv[2])
maxlen = int(sys.argv[3]) if len(sys.argv) > 3 else None
rs = RootSystem(T, N); lam = rs.fundamental(k)
print(f"{typ} om{k}  nakayama nu({k}) = {nakayama(rs)[k]}", flush=True)
nv = nr = bad = 0
for mu, word in sorted(orbit_with_words(rs, lam, max_len=maxlen).items(),
                       key=lambda x: (len(x[1]), x[1])):
    if not word: continue
    v = tuple(v_from_trace(rs, trace(rs, word, lam)))
    X, rows, dv = soc_sigma(rs, k, word)
    r = solve(rs, k, v)
    if not r:
        print(f"  dfs3 no solution for v={v} (word {word})", flush=True); continue
    drows = list(r.values())[0]
    M = max(rows) if rows else 0
    okv = (dv == v); okr = ({M - d: x for d, x in rows.items()} == drows)
    nv += okv; nr += okr
    if not (okv and okr):
        bad += 1
        if bad <= 5:
            print(f"  MISMATCH w={word} v={v} dv={dv}\n     soc={dict(sorted(rows.items()))}"
                  f"\n     dfs3={dict(sorted(drows.items()))}", flush=True)
print(f"  dimvec agree {nv}, rows agree (after reversal) {nr}, mismatches {bad}")
