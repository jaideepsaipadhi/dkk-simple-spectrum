"""Hostile check of the statements added in section sec:nodeone of typeD.tex.

Each test is the statement AS WRITTEN in the paper, not a weaker proxy.
usage: python3 verify_nodeone.py [MAXN=7]
"""
import sys
from rootsys import RootSystem, orbit_with_words
from wordtrace import trace
from weylprops import reduced_words
from socle import soc_sigma

FAST = "--fast" in sys.argv
MAXN = int(sys.argv[1]) if len(sys.argv) > 1 and sys.argv[1].isdigit() else 7
fail = 0
for N in range(4, MAXN + 1):
    rs = RootSystem('D', N); t = N - 2
    n_vs = n_lone = n_cor = n_thm = 0
    for k in rs.I:
        lam = rs.fundamental(k)
        for mu, w in orbit_with_words(rs, lam).items():
            if not w: continue
            tr = trace(rs, w, lam); cs = [(i, c) for _, i, c, _ in tr]
            _, rows, dv = soc_sigma(rs, k, w)
            simple = all(x <= 1 for y in rows.values() for x in y.values())

            # eq:vsumc -- v_i = sum of c_p over occurrences of i
            got = {i: 0 for i in rs.I}
            for i, c in cs: got[i] += c
            if any(got[i] != dv[i-1] for i in rs.I): n_vs += 1

            # lem:loneT -- claimed for 1 <= k <= n-2, over EVERY reduced word
            if k <= N - 2 and not FAST:
                for u in reduced_words(rs, w):
                    if u.count(t) != 1: continue
                    tru = trace(rs, u, lam)
                    cu = [c for _, i, c, _ in tru if i == t]
                    if cu != [1]: n_lone += 1; break

            # cor:nodeone -- ANY fundamental lambda, spin included
            if any(i == t and c >= 2 for i, c in cs) and simple: n_cor += 1

            # thm:typeD in the new form, in its stated range
            if 2 <= k <= N - 2 and simple != (dv[t-1] <= 2): n_thm += 1
    bad = n_vs + n_lone + n_cor + n_thm
    fail += bad
    print(f"D{N}: eq:vsumc {n_vs}, lem:loneT {n_lone}, cor:nodeone {n_cor}, "
          f"thm:typeD(v_t<=2) {n_thm}   -> {'OK' if not bad else 'FAIL'}")
print(("all statements hold" if not fail else f"{fail} FAILURES") +
      f", D4-D{MAXN}, every fundamental lambda, full Weyl orbits")
