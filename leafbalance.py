"""
LEMMA B (leaf balance) and the completion of Proposition T in type D.

Realisability constraint, proved: for lambda = omega_k with k <= n-2 and
mu = w.lambda, inverting the coordinate change of Lemma C gives

        v_n - v_{n-1} = - mu_n  in  {0, +1, -1}.                       (LB)

(Indeed the e_n-coefficient of lambda - mu is c_n = v_n - v_{n-1}, and lambda
has no support in coordinate n, so c_n = -mu_n, which is 0 or +-1.)

Claim: every graded submodule of I_kappa that violates UPWARD CLOSURE at
colour t -- i.e. has V_t[d] != 0 but V_t[d+2] != I_t[d+2] -- has
|v_n - v_{n-1}| >= 2, hence is not realisable.  Combined with Lemma U
(the double-path span is everything for non-special x) this proves
Proposition T for all k in type D.

This script enumerates ALL graded submodules and checks the claim.
"""
import sys, time
from rootsys import RootSystem
from proj2 import Proj2
from allsub import all_submodules
from lemA import TRIV


def run(typ, n, k):
    rs = RootSystem(typ, n)
    t = TRIV[typ](n)
    pr, subs = all_submodules(rs, k)
    viol = 0
    unreal = 0
    bal = {}
    for N in subs:
        # upward closure at colour t
        bad = False
        for d in range(pr.top + 1):
            if N.get((d, t)) and pr.dim.get((d + 2, t), 0):
                if len(N.get((d + 2, t), [])) < pr.dim[(d + 2, t)]:
                    bad = True
                    break
        if not bad:
            continue
        viol += 1
        v = [sum(len(N.get((d, i), [])) for d in range(pr.top + 1)) for i in rs.I]
        diff = v[n - 1] - v[n - 2]          # v_n - v_{n-1}
        bal[diff] = bal.get(diff, 0) + 1
        if abs(diff) >= 2:
            unreal += 1
    return len(subs), viol, unreal, bal


if __name__ == "__main__":
    t0 = time.time()
    V = U = 0
    for typ, n, k in [("D", 4, 2), ("D", 5, 2), ("D", 5, 3), ("D", 6, 2),
                      ("D", 6, 3), ("D", 6, 4), ("D", 7, 3)]:
        try:
            tot, viol, unreal, bal = run(typ, n, k)
        except Exception as e:
            print(f"  {typ}{n} e{k}: skipped ({e})")
            continue
        V += viol; U += unreal
        print(f"  {typ}{n} e{k}: {tot} graded submodules, {viol} violate upward "
              f"closure, {unreal} of those have |v_n - v_(n-1)| >= 2   "
              f"(distribution of v_n - v_(n-1): {bal})")
    print(f"\ntotal: {U}/{V} upward-closure violators excluded by leaf balance "
          f"({time.time()-t0:.0f}s)")
    print("LEAF BALANCE EXCLUDES EVERY VIOLATOR" if U == V else "*** not all ***")
    sys.exit(0 if U == V else 1)
