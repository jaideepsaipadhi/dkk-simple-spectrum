"""
LEMMA I IN TYPE D, PROVED -- verification of the two steps of the argument.

  Lemma I (type D).  Let t = n-2, let i < t be a chain vertex, and let d be a
  degree with dim (I_kappa)_i[d] >= 2.  Then the inward arrow
  (I_kappa)_i[d] -> (I_kappa)_{i+1}[d+1] is injective.

  Proof.  Suppose x != 0 lies in the kernel.  Vertex i has exactly two
  neighbours, i-1 and i+1 (for i < t).  Since I_kappa has SIMPLE SOCLE S_kappa,
  concentrated in the top degree, an element killed by every arrow out of i is
  either zero or that socle generator; as dim (I_kappa)_i[d] >= 2 we are below
  the top degree, so x is not killed by both arrows and

        y := a_{i -> i-1}(x)  !=  0.

  The preprojective relation at i reads
        eps  a_{i+1 -> i} a_{i -> i+1}  +  eps' a_{i-1 -> i} a_{i -> i-1} = 0
  as a map (I_kappa)_i[d] -> (I_kappa)_i[d+2].  Applying it to x and using
  a_{i -> i+1}(x) = 0 gives  a_{i-1 -> i}(y) = 0:  the SAME situation one vertex
  further from t, one degree higher.

  Iterating, we reach a nonzero z in (I_kappa)_1[d + i - 1] killed by the arrow
  1 -> 2.  Vertex 1 of D_n has 2 as its ONLY neighbour, so z is killed by every
  arrow out of 1, hence z lies in soc(I_kappa) = S_kappa.  That forces
  kappa = 1.  But by Lemma P the number of degrees at which a chain vertex i
  carries multiplicity 2 is max(0, i + k - n + 1), which for k = 1 and
  i <= n-2 is 0 -- so dim (I_kappa)_i[d] >= 2 is impossible when kappa = 1.
  Contradiction.  []

This script verifies the two inputs and the conclusion:
  (A) soc(I_kappa) is simple, concentrated in the top degree, at vertex kappa;
  (B) for kappa = 1 no chain vertex carries multiplicity 2;
  (C) the algebraic input to the propagation: at a chain vertex i the two
      double paths through i-1 and through i+1 cancel, which is exactly what
      turns "x killed inward" into "a_{i->i-1}(x) killed inward".

Note the propagation itself has no instances to test: it is an argument by
contradiction, and since the lemma is true the hypothetical kernel elements do
not exist.  What is checkable is (A), (B), (C), and the conclusion, which
inwardlemma.py verifies directly.
"""
import sys, itertools
from rootsys import RootSystem
from proj2 import Proj2
from qvar import rref, P


def socle(pr, rs):
    """{(d,i): dim} of the elements killed by every arrow."""
    out = {}
    for (d, i), m in sorted(pr.dim.items()):
        if m == 0:
            continue
        eqs = []
        for j in rs.adj[i]:
            nj = pr.dim.get((d + 1, j), 0)
            if nj == 0 or (d, i, j) not in pr.arrow:
                continue
            imgs = [pr.act(d, i, [1 if x == b else 0 for x in range(m)], j)
                    for b in range(m)]
            for x in range(nj):
                row = [imgs[b][x] % P for b in range(m)]
                if any(row):
                    eqs.append(row)
        if not eqs:
            out[(d, i)] = m
            continue
        R, piv = rref(eqs, m)
        k = m - len(R)
        if k:
            out[(d, i)] = k
    return out


def check(n):
    rs = RootSystem("D", n)
    t = n - 2
    resA = resB = resC = True
    nA = nB = nC = 0
    for k in rs.I:
        pr = Proj2(rs, k, dmax=60)
        # (A) simple socle in the top degree
        soc = socle(pr, rs)
        nA += 1
        if not (len(soc) == 1 and sum(soc.values()) == 1
                and list(soc)[0][0] == pr.top):
            resA = False
            print(f"    D{n} e{k}: socle = {soc}  *** not simple/top ***")
        # (B) kappa with k = 1: no chain vertex carries multiplicity 2
        if k == 1:
            for (d, i), m in pr.dim.items():
                if i <= t and m >= 2:
                    resB = False
                    print(f"    D{n} e1: vertex {i} degree {d} has mult {m}"
                          f"  *** multiplicity 2 with k=1 ***")
            nB += 1
        # (C) the algebraic input to the propagation: the preprojective
        # relation at a chain vertex i, i.e. that the two double paths
        # i -> i+1 -> i  and  i -> i-1 -> i  are negatives of each other.
        for (d, i), m in sorted(pr.dim.items()):
            if i >= t or m == 0 or i == 1:
                continue
            for b in range(m):
                x = [1 if u == b else 0 for u in range(m)]
                tot = None
                for j in (i - 1, i + 1):
                    if (d, i, j) not in pr.arrow or (d + 1, j, i) not in pr.arrow:
                        continue
                    y = pr.act(d, i, x, j)
                    z = pr.act(d + 1, j, y, i)
                    eps = 1 if i < j else -1
                    z = [(eps * v) % P for v in z]
                    tot = z if tot is None else [(a + c) % P
                                                 for a, c in zip(tot, z)]
                nC += 1
                if tot is not None and any(v % P for v in tot):
                    resC = False
                    print(f"    D{n} e{k}: relation fails at ({d},{i})")
    return resA, resB, resC, nA, nB, nC


if __name__ == "__main__":
    ok = True
    A = B = C = 0
    print("Verifying the ingredients of the type-D proof of Lemma I\n")
    for n in range(4, 10):
        a, b, c, na, nb, nc = check(n)
        ok &= a and b and c
        A += na; B += nb; C += nc
        print(f"  D{n}: (A) simple socle {na} injectives, "
              f"(B) k=1 multiplicity-free on the chain, "
              f"(C) {nc} relation instances   "
              f"{'OK' if a and b and c else 'FAIL'}")
    print(f"\ntotals: {A} injectives, {B} weights with k=1, "
          f"{C} relation instances")
    print("EVERY INGREDIENT VERIFIED" if ok else "*** something failed ***")
    sys.exit(0 if ok else 1)
