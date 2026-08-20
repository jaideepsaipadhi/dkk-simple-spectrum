"""
Toward a PROOF of Proposition T (top truncation) for all k.

Since v = lambda - w.lambda is w-extremal, Savage-Tingley uniqueness says the
graded submodule N of I_kappa with dim N = v is exactly V = soc_sigma(I_kappa).
So Proposition T is a statement about V, and we may use the construction.

Two candidate sandwich bounds, both about DEGREE truncations of I = I_kappa
(recall rad^m(I) = the part in degrees >= m, since Pi is generated in degree
<= 1):

  (U)  V  is contained in  soc^{L}(I),  L = number of blocks = ell(w).
       [each socle step adds at most one socle layer]

  (D)  V  contains  rad^{m}(I)  for m = ell(w_0^lambda) - ell(w).
       [dual statement: I/V is the dual of the extremal module for the
        complementary Weyl element]

This script measures:
  * whether soc^r(I) equals the degree truncation "degrees >= d_max - r + 1"
    (i.e. whether I is a rigid / Loewy module);
  * the minimal m with rad^m(I) contained in V, against ell(w_0^lambda)-ell(w);
  * the maximal r with V contained in soc^r(I), against ell(w);
  * whether (D) alone already forces the colour-t top truncation.
"""
import sys, time
from rootsys import RootSystem, orbit_with_words
from wordtrace import trace
from proj2 import Proj2
from qvar import rref, P
from homological import setup
from lemA import build, TRIV


def rad_in_V(pr, rs, X, m):
    """Is every basis vector of I in degrees >= m already inside X?"""
    for (d, i), dim in pr.dim.items():
        if d >= m and len(X.get((d, i), [])) < dim:
            return False
    return True


def min_rad(pr, rs, X):
    for m in range(pr.top + 2):
        if rad_in_V(pr, rs, X, m):
            return m
    return None


def V_lowest_degree(pr, rs, X):
    ds = [d for (d, i), rows in X.items() if rows]
    return min(ds) if ds else pr.top + 1


def longest_len(rs, k):
    """ell(w_0^lambda) = max number of c-trace steps over the orbit."""
    lam = rs.fundamental(k)
    return max(len(w) for w in orbit_with_words(rs, lam).values())


def run(typ, n, k, cap=300):
    rs, pr, kk = setup(typ, n, k)
    lam = rs.fundamental(k)
    orb = orbit_with_words(rs, lam)
    if len(orb) > cap:
        return None
    L0 = longest_len(rs, k)
    tot = 0
    okD = okU = 0
    for mu, word in orb.items():
        if not word:
            continue
        X = build(rs, pr, word, lam)
        tot += 1
        m = min_rad(pr, rs, X)
        okD += (m is not None and m <= L0 - len(word) + pr.top - L0 + 1 or
                m == pr.top - len(word) + 1)
        lo = V_lowest_degree(pr, rs, X)
        okU += (lo >= pr.top - len(word) + 1)
    return tot, okD, okU, L0, pr.top


if __name__ == "__main__":
    t0 = time.time()
    print("case        elts   V inside deg >= top-l(w)+1    rad^{top-l(w)+1} inside V"
          "      top  l(w_0^lam)")
    T = U = D = 0
    for typ, n, k in [("D", 4, 2), ("D", 5, 2), ("D", 5, 3), ("D", 6, 3),
                      ("D", 6, 4), ("D", 7, 3), ("E", 6, 2), ("E", 6, 3),
                      ("E", 7, 1), ("A", 5, 3)]:
        r = run(typ, n, k)
        if r is None:
            continue
        tot, okD, okU, L0, top = r
        T += tot; U += okU; D += okD
        print(f"{typ}{n} om{k:<3} {tot:5d}   {okU:5d}/{tot:<5d}                "
              f"{okD:5d}/{tot:<5d}                  {top:3d}  {L0:3d}")
    print(f"\ntotals  upper {U}/{T}   lower {D}/{T}   ({time.time()-t0:.0f}s)")
