"""
THE CLASSIFICATION RUN.

For each w in the orbit W.lambda: compute v from the c-trace, solve for the
graded quotient of Pi e_k with that dimension vector, and test SIMPLE SPECTRUM
(every graded piece of dimension <= 1).  Then correlate with combinatorial
invariants of w to find the criterion.
"""
import sys, time
from rootsys import RootSystem, orbit_with_words
from wordtrace import trace, v_from_trace
from weylprops import is_fully_commutative
from qsearch2 import search2

typ = sys.argv[1] if len(sys.argv) > 3 else "D"
n = int(sys.argv[2]) if len(sys.argv) > 3 else 5
k = int(sys.argv[3]) if len(sys.argv) > 3 else 2
dmax = int(sys.argv[4]) if len(sys.argv) > 4 else 8

rs = RootSystem(typ, n); lam = rs.fundamental(k)
orb = orbit_with_words(rs, lam)
print(f"{rs.name()} lambda=omega_{k}   |orbit| = {len(orb)}\n")

rows = []
t0 = time.time()
for mu, word in sorted(orb.items(), key=lambda t: len(t[1])):
    if not word: continue
    v = v_from_trace(rs, trace(rs, word, lam))
    cs = [c for _, _, c, _ in trace(rs, word, lam)]
    try:
        sols = search2(rs, k, v, dmax=dmax)
    except Exception as e:
        print(f"  w={''.join('s'+str(i) for i in word)}  ERROR {e}"); continue
    if not sols:
        status = "UNRESOLVED"
        ss = None
    else:
        gr = list(sols.values())[0]
        ss = all(c <= 1 for r in gr.values() for c in r.values())
        status = "simple" if ss else "NOT simple"
    fc = is_fully_commutative(rs, word) if len(word) <= 11 else None
    rows.append((word, tuple(v), max(cs), fc, ss, status, len(sols)))
    print(f"  w={''.join('s'+str(i) for i in word):<26} v={tuple(v)}  "
          f"cmax={max(cs)}  FC={str(fc):<5} -> {status}"
          + ("" if len(sols) <= 1 else f"  ({len(sols)} solutions!)"))

print(f"\nelapsed {time.time()-t0:.0f}s")
good = [r for r in rows if r[4] is True]
bad  = [r for r in rows if r[4] is False]
unk  = [r for r in rows if r[4] is None]
print(f"simple spectrum: {len(good)}   NOT simple: {len(bad)}   unresolved: {len(unk)}")
print("\n--- correlation with c_max ---")
for cm in sorted({r[2] for r in rows}):
    g = sum(1 for r in rows if r[2]==cm and r[4] is True)
    b = sum(1 for r in rows if r[2]==cm and r[4] is False)
    print(f"  c_max={cm}: simple={g}  not simple={b}")
print("\n--- correlation with fully commutative ---")
for f in (True, False):
    g = sum(1 for r in rows if r[3] is f and r[4] is True)
    b = sum(1 for r in rows if r[3] is f and r[4] is False)
    print(f"  FC={f}: simple={g}  not simple={b}")
if bad:
    print("\nNOT simple examples:")
    for r in bad[:10]:
        print(f"   w={''.join('s'+str(i) for i in r[0])}  v={r[1]}  cmax={r[2]}  FC={r[3]}")
