"""Among the E6 lam=om4 elements with a c-trace entry 3, what separates the unique
simple-spectrum one from the rest? The verdict depends only on v, so group by v first."""
import collections
from rootsys import RootSystem, orbit_with_words
from wordtrace import trace, v_from_trace
from dfs3 import solve
rs = RootSystem("E", 6); k = 4; lam = rs.fundamental(k)
orb = orbit_with_words(rs, lam)
rec = []
for mu, word in orb.items():
    if not word: continue
    tr = trace(rs, word, lam)
    cs = [c for _, _, c, _ in tr]
    if max(cs) != 3: continue
    rec.append((tuple(v_from_trace(rs, tr)), cs, len(word), tuple(mu)))
vs = sorted({r[0] for r in rec})
print(f"{len(rec)} elements with maxc=3, {len(vs)} distinct v")
verdict = {}
for v in vs:
    r = solve(rs, k, v)
    if not r: verdict[v] = None; continue
    rows = list(r.values())[0]
    verdict[v] = all(c <= 1 for x in rows.values() for c in x.values())
rec = [r for r in rec if verdict[r[0]] is not None]
print(f"{len(rec)} with a solution, {sum(verdict[r[0]] for r in rec)} simple")
feats = lambda cs, L: dict(n3=cs.count(3), n2=cs.count(2), last3=cs[-1] == 3,
                           only_last=(cs.count(3) == 1 and cs[-1] == 3 and cs.count(2) == 0), L=L)
for feat in ("n3", "n2", "last3", "only_last", "L"):
    t = collections.Counter(); s = collections.Counter()
    for v, cs, L, mu in rec:
        k_ = feats(cs, L)[feat]; t[k_] += 1; s[k_] += verdict[v]
    print(f"\n{feat}: value -> (total, simple)")
    for val in sorted(t, key=str):
        print(f"   {val!s:7} -> {t[val]:5} {s[val]:3}"
              + ("  <== all simple" if s[val] == t[val] else ""))
print("\nthe simple ones:")
for v, cs, L, mu in rec:
    if verdict[v]: print(f"   L={L} mu={mu} v={v} c={cs}")
