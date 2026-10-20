"""Does the c-trace alone decide simple spectrum?  Tested over all of E6 and E7.
Two candidate statements:
  (A) maxc >= 4  ==>  never simple spectrum.
  (B) maxc == 3  ==>  simple spectrum  <=>  the c-trace is (1,...,1,3).
usage: python3 echar.py E6 [E7 ...]"""
import sys, collections, pickle, os
from rootsys import RootSystem, orbit_with_words
from wordtrace import trace, v_from_trace
from socle import soc_sigma

def records(typ):
    import os as _os; cache = f".echar_{typ}_{_os.environ.get(chr(68)+chr(75)+chr(75)+chr(80),7)}.pkl"
    if os.path.exists(cache): return pickle.load(open(cache, "rb"))
    T, N = typ[0], int(typ[1:]); rs = RootSystem(T, N); out = []
    for k in rs.I:
        lam = rs.fundamental(k)
        for mu, word in orbit_with_words(rs, lam).items():
            if not word: continue
            cs = [c for _, _, c, _ in trace(rs, word, lam)]
            X, rows, dv = soc_sigma(rs, k, word)
            ss = all(c <= 1 for x in rows.values() for c in x.values())
            out.append((k, tuple(cs), ss, dv))
        print(f"  {typ} om{k} done", flush=True)
    pickle.dump(out, open(cache, "wb"))
    return out

rec = []
for typ in sys.argv[1:]: rec += [(typ,) + r for r in records(typ)]
print(f"\n{len(rec)} elements total")
print("\nmaxc -> (count, simple count)")
t = collections.Counter(); s = collections.Counter()
for typ, k, cs, ss, dv in rec:
    m = max(cs); t[m] += 1; s[m] += ss
for m in sorted(t): print(f"   maxc={m}: {t[m]:7} {s[m]:5}")
print("\n(A) any element with maxc >= 4 and simple spectrum?")
bad = [r for r in rec if max(r[2]) >= 4 and r[3]]
print(f"    {len(bad)} found" + ("" if not bad else f"  e.g. {bad[0]}"))
print("\n(B) restricted to maxc == 3:  simple  vs  trace == (1,...,1,3)")
tab = collections.Counter()
for typ, k, cs, ss, dv in rec:
    if max(cs) != 3: continue
    cond = (cs[-1] == 3 and all(c == 1 for c in cs[:-1]))
    tab[(cond, ss)] += 1
for key in sorted(tab): print(f"    trace=(1..1,3): {key[0]!s:5}  simple: {key[1]!s:5} -> {tab[key]}")
print("\n(B') same test for maxc == 2:")
tab = collections.Counter()
for typ, k, cs, ss, dv in rec:
    if max(cs) != 2: continue
    cond = (cs[-1] == 2 and all(c == 1 for c in cs[:-1]))
    tab[(cond, ss)] += 1
for key in sorted(tab): print(f"    trace=(1..1,2): {key[0]!s:5}  simple: {key[1]!s:5} -> {tab[key]}")
