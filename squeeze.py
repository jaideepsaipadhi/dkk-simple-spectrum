"""
THE SQUEEZE:  dim I_kappa = lambda - w_0 lambda  (root coordinates),
lambda = omega_kappa, I_kappa = injective envelope of S_kappa over Pi.

This is what upgrades the Crawley-Boevey identity

        dim soc_i(I/X) = <alpha_i^vee, lambda - dim X> + dim top_i(X)   (CB)

into a PROOF of Lemma S.  Extend a reduced word for w to a reduced word for
the longest element w_0 and run the c-trace.  By (CB) each step raises dim X
by at least its c, so

        dim soc_{sigma_full}(I) >= sum of all c = lambda - w_0 lambda.

But soc_{sigma_full}(I) <= I and dim I = lambda - w_0 lambda, so equality
holds and EVERY step contributes exactly c: all the top-terms vanish.  In
particular this holds along the prefix corresponding to w, which is Lemma S.

This script checks the one input to that argument.
"""
import sys
from rootsys import RootSystem, orbit_with_words
from proj2 import Proj2


def socv(rs, j):
    pr = Proj2(rs, j)
    return max((i for (d, i), m in pr.dim.items() if d == pr.top and m), default=None)


def longest_weight(rs, k):
    """w_0 . omega_k = the unique antidominant element of the orbit."""
    lam = rs.fundamental(k)
    orb = orbit_with_words(rs, lam)
    for mu in orb:
        if all(x <= 0 for x in mu):
            return mu
    raise RuntimeError("no antidominant weight found")


def check(typ, n):
    rs = RootSystem(typ, n)
    nu = {j: socv(rs, j) for j in rs.I}
    inv = {b: a for a, b in nu.items()}
    print(f"=== {rs.name()}   Nakayama nu = {nu}")
    ok = True
    for k in rs.I:
        pr = Proj2(rs, inv[k], dmax=60)
        d = tuple(sum(pr.dim.get((dd, i), 0) for dd in range(pr.top + 1))
                  for i in rs.I)
        lam = rs.fundamental(k)
        w0lam = longest_weight(rs, k)
        target = tuple(int(x) for x in
                       rs.root_coords([lam[j - 1] - w0lam[j - 1] for j in rs.I]))
        good = d == target
        ok &= good
        print(f"   k={k}: dim I_k = {d}   lambda - w0.lambda = {target}   "
              f"{'OK' if good else 'FAIL'}")
    return ok


if __name__ == "__main__":
    allok = True
    for typ, n in [("A", 3), ("A", 4), ("A", 5), ("D", 4), ("D", 5), ("D", 6),
                   ("D", 7), ("E", 6), ("E", 7)]:
        allok &= check(typ, n)
    print("\n" + ("dim I_kappa = lambda - w_0 lambda IN EVERY CASE"
                  if allok else "*** FAILS ***"))
    sys.exit(0 if allok else 1)
