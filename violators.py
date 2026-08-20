"""
Structure of the upward-closure violators.

Goal: prove that a graded submodule N of I_kappa whose colour-t profile is not
a top truncation has delta(v) = 2 v_kappa - v^T C v > 0, i.e. that its
dimension vector is not w-extremal.

The right framework is the MV polytope.  Baumann-Kamnitzer (Represent. Theory
16 (2012)) show that Pol(T) is the convex hull of the dimension vectors of the
submodules of T, and Baumann-Kamnitzer-Tingley (Publ. IHES 120 (2014), Sect.
3.1) prove RIGIDITY: if x is a vertex of Pol(T) then T has a UNIQUE subobject
with dimension vector x.  For T = I_kappa the vertices are exactly the
w-extremal dimension vectors.  Hence:

    two distinct submodules with the same dimension vector v
        =>  v is not a vertex  =>  v is not extremal  =>  delta(v) > 0.

So it suffices to exhibit a SECOND submodule with the same dimension vector.
The same framework explains the numerology: since dim U lies in Pol(I_kappa),
mu = lambda - dim U is a weight of L(omega_k), so in type D it has entries in
{0,+-1} with |supp mu| congruent to k mod 2, and delta = k - |supp mu|.

Two candidate constructions for the second submodule:

  (tau) the diagram automorphism swapping the two leaves n-1, n.  It fixes
        I_kappa when kappa <= n-2, and sends a submodule N to another submodule
        with the leaf coordinates of v exchanged.  So when v_{n-1} = v_n it
        preserves the dimension vector, and it moves N whenever N involves one
        of the two special lines at t asymmetrically.

  (any) brute force: is there any other graded submodule with the same v?

This script measures, for every violator: v, the corresponding mu in
epsilon-coordinates, |supp mu| against k, delta, whether v_{n-1} = v_n, whether
tau moves N, and the total number of graded submodules sharing v.
"""
import sys
from collections import Counter, defaultdict
from rootsys import RootSystem
from proj2 import Proj2
from qvar import rref, P
from allsub import all_submodules
from epsilon import eps_of_mu
from lemA import TRIV


def violates(pr, N, t):
    for d in range(pr.top + 1):
        if N.get((d, t)) and pr.dim.get((d + 2, t), 0):
            if len(N.get((d + 2, t), [])) < pr.dim[(d + 2, t)]:
                return True
    return False


def dimvec(pr, rs, N):
    return tuple(sum(len(N.get((d, i), [])) for d in range(pr.top + 1))
                 for i in rs.I)


def delta(rs, k, v):
    return 2 * v[k - 1] - sum(v[i - 1] * rs.A[i - 1][j - 1] * v[j - 1]
                              for i in rs.I for j in rs.I)


def run(typ, n, k):
    rs = RootSystem(typ, n)
    t = TRIV[typ](n)
    pr, subs = all_submodules(rs, k)
    bydim = defaultdict(list)
    for N in subs:
        bydim[dimvec(pr, rs, N)].append(N)
    rows = []
    for N in subs:
        if not violates(pr, N, t):
            continue
        v = dimvec(pr, rs, N)
        e = eps_of_mu(n, k, list(v))
        supp = sum(1 for x in e if x)
        rows.append(dict(v=v, eps=tuple(e), supp=supp, k=k,
                         delta=delta(rs, k, v),
                         leafeq=(v[n - 2] == v[n - 1]),
                         nsub=len(bydim[v])))
    return rows


if __name__ == "__main__":
    CASES = [("D", 5, 3), ("D", 6, 3), ("D", 6, 4), ("D", 7, 3)]
    allrows = []
    for typ, n, k in CASES:
        rows = run(typ, n, k)
        allrows += rows
        c = Counter((r["delta"], r["supp"], r["leafeq"], r["nsub"] > 1)
                    for r in rows)
        print(f"  {typ}{n} e{k}: {len(rows)} violators")
        for key, cnt in sorted(c.items()):
            d, sp, le, multi = key
            print(f"      delta={d}  |supp mu|={sp} (k={k})  "
                  f"v_(n-1)=v_n: {le}  another submodule with same v: {multi}"
                  f"   x{cnt}")
    print()
    print(f"total violators: {len(allrows)}")
    print(f"  all have delta > 0                : "
          f"{all(r['delta'] > 0 for r in allrows)}")
    print(f"  all have |supp mu| < k            : "
          f"{all(r['supp'] < r['k'] for r in allrows)}")
    print(f"  all have delta = k - |supp mu|    : "
          f"{all(r['delta'] == r['k'] - r['supp'] for r in allrows)}")
    print(f"  all share v with another submodule: "
          f"{all(r['nsub'] > 1 for r in allrows)}")
    print(f"  fraction with v_(n-1) = v_n       : "
          f"{sum(1 for r in allrows if r['leafeq'])}/{len(allrows)}")
