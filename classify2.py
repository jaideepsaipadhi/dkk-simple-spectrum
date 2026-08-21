import sys, time
from rootsys import RootSystem, orbit_with_words
from wordtrace import trace, v_from_trace
from weylprops import is_fully_commutative
from dfs import all_quotient_dims

typ = sys.argv[1] if len(sys.argv) > 3 else "D"
n = int(sys.argv[2]) if len(sys.argv) > 3 else 5
k = int(sys.argv[3]) if len(sys.argv) > 3 else 2
dmax = int(sys.argv[4]) if len(sys.argv) > 4 else 8
rs = RootSystem(typ, n); lam = rs.fundamental(k)
tri = [i for i in rs.I if len(rs.adj[i]) >= 3]
orb = orbit_with_words(rs, lam)
print(f"{rs.name()} om{k}  |orbit|={len(orb)}  trivalent={tri}")
rows = []; t0 = time.time()
for mu, word in sorted(orb.items(), key=lambda t: len(t[1])):
    if not word: continue
    v = v_from_trace(rs, trace(rs, word, lam))
    cs = [c for _, _, c, _ in trace(rs, word, lam)]
    res = all_quotient_dims(rs, k, v, dmax=dmax)
    if not res:
        rows.append((word, v, max(cs), None, None, 0)); continue
    profs = list(res.values())
    ss = all(all(c <= 1 for r in p.values() for c in r.values()) for p in profs)
    fc = is_fully_commutative(rs, word) if len(word) <= 12 else None
    rows.append((word, v, max(cs), fc, ss, len(profs)))
print(f"elapsed {time.time()-t0:.0f}s")
un = [r for r in rows if r[4] is None]
print(f"resolved {len(rows)-len(un)}/{len(rows)}   unresolved {len(un)}")
multi = [r for r in rows if r[5] > 1]
print(f"elements with >1 graded profile: {len(multi)}")
g = [r for r in rows if r[4] is True]; b = [r for r in rows if r[4] is False]
print(f"simple {len(g)}   not simple {len(b)}")
print("\nc_max:")
for cm in sorted({r[2] for r in rows if r[4] is not None}):
    print(f"  c_max={cm}: simple={sum(1 for r in g if r[2]==cm)} not={sum(1 for r in b if r[2]==cm)}")
print("FC:")
for f in (True, False):
    print(f"  FC={f}: simple={sum(1 for r in g if r[3] is f)} not={sum(1 for r in b if r[3] is f)}")
print("trivalent v_t values:")
for t in tri:
    print(f"  node {t}: simple v_t in {sorted({r[1][t-1] for r in g})}, "
          f"not-simple v_t in {sorted({r[1][t-1] for r in b})}")
