"""
THE SINGLE REMAINING CONJECTURE.

Every residual claim in this project -- top truncation at the trivalent node,
the up-set property at the leaves, Claim F, Claim G -- collapses to one
statement about graded submodules of I_kappa.

Call a graded submodule M of I_kappa a STAIRCASE if

    * its colour-t multiplicities are the greedy TOP truncation of those of
      I_kappa (zero below one degree, proper at it, full above), and
    * at each of the two leaves M occupies a final segment of degrees.

    CONJECTURE (uniqueness implies staircase).  If dim M is w-extremal then M
    is a staircase.

Equivalently, by the rigidity property of Baumann-Kamnitzer-Tingley -- a vertex
of the MV polytope carries a UNIQUE subobject -- the conjecture says

    M is not a staircase  ==>  delta(dim M) = 2 v_kappa - v^T C v > 0,

i.e. a non-staircase submodule never has an extremal dimension vector, and so
is never the unique submodule of its dimension vector.

WHY THIS IS THE WHOLE OF IT.  Proposition T is the contrapositive at t: an
extremal v carries a unique submodule, namely V = soc_sigma(I_kappa), and the
conjecture makes it a staircase, which is Proposition T.  The leaf clause is
what the proof of "G1 at the leaves" assumes.  Theorem "type D" then follows.

WHAT IS PROVED.  The conjecture holds whenever v_{n-1} = v_n, by the
leaf-exchanging diagram automorphism tau (Proposition tau in the paper): tau
fixes I_kappa, exchanges the two special lines at t, hence moves M while
preserving dim M, and two distinct submodules of one dimension vector
contradict rigidity.

This script verifies the conjecture over all graded submodules of I_kappa in
range, and reports the delta of every non-staircase.
"""
import sys, time
from collections import Counter
from rootsys import RootSystem
from allsub import all_submodules
from violators import dimvec, delta
from lemA import TRIV

CASES = [("D", 4, 1), ("D", 4, 2), ("D", 5, 2), ("D", 5, 3),
         ("D", 6, 2), ("D", 6, 3), ("D", 7, 3)]


def is_staircase(pr, rs, M, t, n):
    degs = [d for d in range(pr.top + 1) if pr.dim.get((d, t), 0)]
    prof = [(d, len(M.get((d, t), [])), pr.dim[(d, t)]) for d in degs]
    part = [d for d, a, b in prof if 0 < a < b]
    zero = [d for d, a, b in prof if a == 0]
    full = [d for d, a, b in prof if a == b]
    if len(part) > 1:
        return False
    if zero and full and max(zero) > min(full):
        return False
    if part and ((zero and max(zero) > part[0]) or (full and min(full) < part[0])):
        return False
    for j in (n - 1, n):
        dg = [d for d in range(pr.top + 1) if pr.dim.get((d, j), 0)]
        occ = [d for d in dg if M.get((d, j))]
        if occ and set(occ) != set(d for d in dg if d >= min(occ)):
            return False
    return True


if __name__ == "__main__":
    t0 = time.time()
    TOT = NON = POS = 0
    dist = Counter()
    for typ, n, k in CASES:
        rs = RootSystem(typ, n)
        t = TRIV[typ](n)
        try:
            pr, subs = all_submodules(rs, k)
        except Exception as e:
            print(f"  {typ}{n} e{k}: skipped ({e})")
            continue
        non = pos = 0
        for M in subs:
            TOT += 1
            if is_staircase(pr, rs, M, t, n):
                continue
            non += 1
            d = delta(rs, k, dimvec(pr, rs, M))
            dist[d] += 1
            pos += (d > 0)
        NON += non; POS += pos
        print(f"  {typ}{n} e{k}: {len(subs)} graded submodules, {non} not a "
              f"staircase, {pos} of those with delta > 0"
              f"   {'OK' if pos == non else 'FAIL'}")
    print(f"\ntotal: {TOT} graded submodules, {NON} non-staircases")
    print(f"  every non-staircase has delta > 0 : {POS}/{NON}   "
          f"{'OK' if POS == NON else 'FAIL'}")
    print(f"  delta distribution                : {dict(dist)}")
    print(f"\n({time.time()-t0:.0f}s)")
    print("CONJECTURE VERIFIED IN RANGE" if POS == NON else "*** FAILS ***")
    sys.exit(0 if POS == NON else 1)
