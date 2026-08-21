"""
HYPOTHESIS from the D4/D5 data:

    simple spectrum  <=>  v_t <= 2, where t is the TRIVALENT node.

D4 (t=2): not-simple all have v_2 = 3; simple all have v_2 <= 2.
D5 (t=3): not-simple all have v_3 in {3,4}; simple all have v_3 <= 2.
Note it is specifically the trivalent node: D5 has simple elements with
v_2 = 3, e.g. Ex 5.11 v=(1,3,2,1,1).
"""
import sys
from rootsys import RootSystem, orbit_with_words
from wordtrace import trace, v_from_trace
from qsearch2 import search2

typ = sys.argv[1] if len(sys.argv) > 3 else "D"
n = int(sys.argv[2]) if len(sys.argv) > 3 else 5
k = int(sys.argv[3]) if len(sys.argv) > 3 else 2
dmax = int(sys.argv[4]) if len(sys.argv) > 4 else 8
rs = RootSystem(typ, n); lam = rs.fundamental(k)
tri = [i for i in rs.I if len(rs.adj[i]) >= 3]
print(f"{rs.name()} om{k}   trivalent node(s): {tri}")

agree = dis = unres = 0
bad = []
for mu, word in sorted(orbit_with_words(rs, lam).items(), key=lambda t: len(t[1])):
    if not word: continue
    v = v_from_trace(rs, trace(rs, word, lam))
    try:
        sols = search2(rs, k, v, dmax=dmax)
    except Exception:
        unres += 1; continue
    if not sols:
        unres += 1; continue
    gr = list(sols.values())[0]
    truth = all(c <= 1 for r in gr.values() for c in r.values())
    pred = all(v[t-1] <= 2 for t in tri)
    if truth == pred: agree += 1
    else:
        dis += 1
        bad.append((word, v, truth, pred))
print(f"  hypothesis agrees {agree}, disagrees {dis}, unresolved {unres}")
for w, v, t, p in bad[:8]:
    print(f"    w={''.join('s'+str(i) for i in w)}  v={v}  truth={t} pred={p}")
