"""
LEMMA I (inward injectivity on the injective module) -- this proves
"it is enough to check the trivalent node", uniformly and in particular in
type E, where it is an EXHAUSTIVE finite computation.

  Lemma I.  Let t be the trivalent node and, for i != t, let nx(i) be the
  neighbour of i one step closer to t.  Then for every fundamental weight
  kappa and every degree d with  dim (I_kappa)_i[d] >= 2,  the arrow

        (I_kappa)_i[d]  -->  (I_kappa)_{nx(i)}[d+1]

  is injective.

  Corollary.  Let V be any graded submodule of I_kappa (in particular the
  extremal fixed point).  If dim V_i[d] >= 2 for some i != t, then, since V is
  a submodule, the arrow carries V_i[d] into V_{nx(i)}[d+1] and, being
  injective on the ambient (I_kappa)_i[d] of dimension >= 2, it is injective
  on V_i[d].  Hence dim V_{nx(i)}[d+1] >= 2.  Iterating dist(i,t) times,
  dim V_t[d + dist(i,t)] >= 2.  Contrapositive: if V_t is multiplicity-free
  then so is V.

For E6, E7 and E8 there are only finitely many triples (kappa, i, d), so the
verification below is complete rather than a sample.  Types D_n and A_n form
families; type D is covered independently by the window-containment argument,
and in type A the question is vacuous (every fundamental weight is minuscule).
"""
import sys, time
from rootsys import RootSystem
from qvar import rref, P
from proj2 import Proj2

TRIV = {"D": lambda n: n - 2, "E": lambda n: 4}
DMAX = 60


def dist_to(rs, t):
    D, front = {t: 0}, [t]
    while front:
        nf = []
        for x in front:
            for y in rs.adj[x]:
                if y not in D:
                    D[y] = D[x] + 1
                    nf.append(y)
        front = nf
    return D


def run(typ, n):
    rs = RootSystem(typ, n)
    t = TRIV[typ](n)
    D = dist_to(rs, t)
    nx = {i: min(rs.adj[i], key=lambda y: D[y]) for i in rs.I if i != t}
    tot = ok = 0
    bad = []
    for k in rs.I:
        pr = Proj2(rs, k, dmax=DMAX)
        assert pr.top < DMAX, f"module did not terminate for {typ}{n} e{k}"
        for (d, i), m in sorted(pr.dim.items()):
            if i == t or m < 2:
                continue
            j = nx[i]
            nj = pr.dim.get((d + 1, j), 0)
            if nj == 0:
                bad.append((k, i, d, m, 0))
                tot += 1
                continue
            imgs = []
            for b in range(m):
                e = [1 if x == b else 0 for x in range(m)]
                imgs.append(pr.act(d, i, e, j) if (d, i, j) in pr.arrow else [0] * nj)
            R, _ = rref([y[:] for y in imgs], nj)
            tot += 1
            if len(R) == m:
                ok += 1
            else:
                bad.append((k, i, d, m, len(R)))
    return tot, ok, bad


if __name__ == "__main__":
    t0 = time.time()
    allok = True
    print("Type E -- exhaustive: every fundamental weight, vertex and degree\n")
    TE = OE = 0
    for typ, n in [("E", 6), ("E", 7), ("E", 8)]:
        tot, ok, bad = run(typ, n)
        TE += tot; OE += ok
        allok &= (ok == tot)
        print(f"  {typ}{n}: {ok}/{tot} graded pieces of multiplicity >= 2 have "
              f"injective inward arrow   {'OK' if ok == tot else bad[:4]}")
    print(f"\n  type E total: {OE}/{TE}")

    print("\nType D -- same check, for the ranks that fit (D is also covered")
    print("independently by window containment)\n")
    TD = OD = 0
    for n in range(4, 11):
        tot, ok, bad = run("D", n)
        TD += tot; OD += ok
        allok &= (ok == tot)
        print(f"  D{n}: {ok}/{tot}   {'OK' if ok == tot else bad[:4]}")
    print(f"\n  type D total: {OD}/{TD}")

    print(f"\n({time.time()-t0:.0f}s)")
    print("LEMMA I VERIFIED" if allok else "*** LEMMA I FAILS ***")
    sys.exit(0 if allok else 1)
