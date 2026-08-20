"""
LEMMA I (inward injectivity) -- a candidate proof of "checking t suffices".

For a vertex i != t let  nx(i)  be the neighbour of i one step closer to the
trivalent node t.  Candidate:

        the arrow  V_i[d] --> V_{nx(i)}[d+1]  is injective for all d.

If so, iterating along the arm gives  dim V_i[d] <= dim V_t[d + dist(i,t)],
so a multiplicity-free colour-t part forces V to be multiplicity-free
everywhere -- Proposition "check only t", uniformly in the type, with no
window-containment argument and no case analysis.

This script tests it, and if it fails, reports the weaker statement
dim V_i[d] <= dim V_t[d + dist(i,t)] directly.
"""
import sys, time
from collections import Counter
from rootsys import RootSystem, orbit_with_words
from wordtrace import trace
from qvar import rref, P
from proj2 import Proj2
from homological import setup
from lemA import build, TRIV


def dist_to(rs, t):
    D = {t: 0}
    front = [t]
    while front:
        nf = []
        for x in front:
            for y in rs.adj[x]:
                if y not in D:
                    D[y] = D[x] + 1
                    nf.append(y)
        front = nf
    return D


def run(typ, n, k, cap=400):
    rs, pr, kk = setup(typ, n, k)
    t = TRIV[typ](n)
    D = dist_to(rs, t)
    nx = {i: min(rs.adj[i], key=lambda y: D[y]) for i in rs.I if i != t}
    lam = rs.fundamental(k)
    orb = orbit_with_words(rs, lam)
    if len(orb) > cap:
        return None
    inj_tot = inj_ok = 0
    dom_tot = dom_ok = 0
    for mu, word in orb.items():
        if not word:
            continue
        X = build(rs, pr, word, lam)
        for d in range(pr.top + 1):
            for i in rs.I:
                rows = X.get((d, i), [])
                if not rows or i == t:
                    continue
                j = nx[i]
                nj = pr.dim.get((d + 1, j), 0)
                imgs = [pr.act(d, i, r, j) for r in rows] if nj else []
                R, _ = rref([y[:] for y in imgs], nj) if nj else ([], [])
                inj_tot += 1
                inj_ok += (len(R) == len(rows))
                dom_tot += 1
                dom_ok += (len(rows) <= len(X.get((d + D[i], t), [])) if
                           X.get((d + D[i], t)) is not None else len(rows) <= 0)
    return inj_tot, inj_ok, dom_tot, dom_ok


CASES = [("D", 4, 2), ("D", 5, 2), ("D", 5, 3), ("D", 6, 3),
         ("E", 6, 2), ("E", 6, 3), ("E", 6, 5), ("E", 7, 1), ("E", 7, 7)]

if __name__ == "__main__":
    t0 = time.time()
    IT = IO = DT = DO = 0
    for typ, n, k in CASES:
        r = run(typ, n, k)
        if r is None:
            continue
        a, b, c, e = r
        IT += a; IO += b; DT += c; DO += e
        print(f"  {typ}{n} om{k}: inward arrow injective {b}/{a};  "
              f"dim V_i[d] <= dim V_t[d+dist] {e}/{c}")
    print(f"\ntotals: injectivity {IO}/{IT},  domination {DO}/{DT}   "
          f"({time.time()-t0:.0f}s)")
    sys.exit(0 if IO == IT else 1)
