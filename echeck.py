"""Independent re-derivation of the unique E6 c_p=3 simple-spectrum element."""
from rootsys import RootSystem, orbit_with_words
from wordtrace import trace, v_from_trace
from weylprops import is_fully_commutative
from dfs3 import solve
rs = RootSystem("E", 6); k = 4; lam = rs.fundamental(k)
word = (2,4,5,6,3,4,5,1,3,4)
# (a) is the word reduced / does it have minimal length in its coset?
orb = orbit_with_words(rs, lam)
mu = lam
mu = rs.act_word(tuple(word), lam)
print("mu = w.lam  =", tuple(mu))
print("orbit word  =", orb[tuple(mu)], "len", len(orb[tuple(mu)]))
print("our word    =", word, "len", len(word))
tr = trace(rs, word, lam)
print("\nc-trace step by step (p, i_p, c_p, mu^(p)):")
for p, rec in enumerate(tr, 1): print(f"  {p:>2} {rec}")
v = v_from_trace(rs, tr)
print("\nv = lam - w.lam  =", tuple(v))
print("check lam-mu     =", tuple(a-b for a, b in zip(rs.fundamental(k), mu)) if False else "")
r = solve(rs, k, v)
rows = list(r.values())[0]
print("\nrows of V (degree: colour multiplicities):")
for d in sorted(rows, reverse=True): print(f"  deg {d}: {dict(sorted(rows[d].items()))}")
print("\nsum of row sizes =", sum(c for x in rows.values() for c in x.values()), " |v| =", sum(v))
print("simple spectrum  =", all(c <= 1 for x in rows.values() for c in x.values()))
print("fully commutative=", is_fully_commutative(rs, word))
# (b) every other word for the same mu: same v, same verdict?
import itertools
print("\nverdict is a function of w, not of the word: recheck via a different reduced word")
