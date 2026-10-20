"""The complete list of type-E extremal weights with a c-trace entry 3 and simple spectrum
(appendix item (3)).  Exhaustive over E6 and E7, all fundamental lambda."""
from rootsys import RootSystem, orbit_with_words
from wordtrace import trace, v_from_trace
from weylprops import is_fully_commutative
from socle import soc_sigma
out = []
for typ in ("E6", "E7"):
    T, N = typ[0], int(typ[1:]); rs = RootSystem(T, N)
    for k in rs.I:
        lam = rs.fundamental(k)
        for mu, word in orbit_with_words(rs, lam).items():
            if not word: continue
            cs = [c for _, _, c, _ in trace(rs, word, lam)]
            if max(cs) != 3: continue
            X, rows, dv = soc_sigma(rs, k, word)
            if not all(c <= 1 for x in rows.values() for c in x.values()): continue
            out.append((typ, k, word, tuple(mu), dv, tuple(cs), rows,
                        is_fully_commutative(rs, word)))
print(f"{len(out)} elements\n")
for n, (typ, k, word, mu, v, cs, rows, fc) in enumerate(
        sorted(out, key=lambda z: (z[0], z[1], len(z[2]))), 1):
    M = max(rows)
    shape = " | ".join("{" + ",".join(str(c) for c in sorted(rows[M-d])
                                      for _ in range(rows[M-d][c])) + "}"
                       for d in range(M, -1, -1) if (M-d) in rows)
    print(f"({n}) {typ}  lam=om{k}  l(w)={len(word)}  fully commutative: {fc}")
    print(f"     w   = {''.join('s'+str(i) for i in word)}")
    print(f"     mu  = {mu}   (omega coordinates)")
    print(f"     v   = {v}   |v| = {sum(v)}")
    print(f"     c   = {list(cs)}")
    print(f"     V   = {shape}   (socle layers, bottom first)")
    print()
