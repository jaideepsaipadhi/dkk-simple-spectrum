"""
THE GREEDY-FROM-THE-TOP MECHANISM.

Key elementary observation: for a graded submodule X of I and a homogeneous
x in I_i[d], one has x in soc_i(I/X) if and only if every arrow i -> j sends x
into X_j[d+1].  In particular:

    if  X_j[d+1] = I_j[d+1]  for every neighbour j of i,
    then  soc_i(I/X)  contains all of  I_i[d].

So socle steps are forced to fill a colour from the top down, as soon as the
neighbouring degrees above have been exhausted.  This script tests the
incremental form of Proposition T:

  (G1)  at every colour-t block of the c-trace, the newly added colour-t
        elements occupy the HIGHEST colour-t degrees of I_kappa that are not
        already full;
  (G2)  hence the running colour-t profile of X is at every stage the greedy
        top truncation of the colour-t profile of I_kappa.

(G2) at the last step is exactly Proposition T.
"""
import sys, time
from rootsys import RootSystem, orbit_with_words
from wordtrace import trace
from proj2 import Proj2
from gls import socle_step
from homological import setup
from lemA import TRIV


def greedy_top(prof, total):
    """prof: list of (degree, mult) increasing; fill from the top."""
    out, left = {}, total
    for d, m in reversed(prof):
        take = min(m, left)
        if take:
            out[d] = take
        left -= take
        if left == 0:
            break
    return out, left


def run(typ, n, k, cap=400):
    rs, pr, kk = setup(typ, n, k)
    t = TRIV[typ](n)
    lam = rs.fundamental(k)
    orb = orbit_with_words(rs, lam)
    if len(orb) > cap:
        return None
    prof = [(d, pr.dim[(d, t)]) for d in range(pr.top + 1) if pr.dim.get((d, t), 0)]
    steps = okstep = 0
    elts = okelt = 0
    for mu, word in orb.items():
        if not word:
            continue
        X = {}
        good = True
        for _, i, c, _ in trace(rs, word, lam):
            for _ in range(c):
                X = socle_step(pr, rs, X, i)
            if i == t:
                got = {d: len(X[(d, t)]) for d in range(pr.top + 1) if X.get((d, t))}
                want, left = greedy_top(prof, sum(got.values()))
                steps += 1
                ok = (got == want and left == 0)
                okstep += ok
                good &= ok
        elts += 1
        okelt += good
    return steps, okstep, elts, okelt, [m for _, m in prof]


CASES = [("D", 4, 2), ("D", 5, 2), ("D", 5, 3), ("D", 6, 2), ("D", 6, 3),
         ("D", 6, 4), ("D", 7, 3), ("E", 6, 2), ("E", 6, 3), ("E", 6, 5),
         ("E", 7, 1), ("E", 7, 7)]

if __name__ == "__main__":
    t0 = time.time()
    S = OS = E = OE = 0
    for typ, n, k in CASES:
        r = run(typ, n, k)
        if r is None:
            print(f"  {typ}{n} om{k}: skipped")
            continue
        s, os_, e, oe, prof = r
        S += s; OS += os_; E += e; OE += oe
        print(f"  {typ}{n} om{k}: (G1/G2) intermediate profiles greedy "
              f"{os_}/{s} colour-t blocks, {oe}/{e} elements   "
              f"(I profile {prof})   {'OK' if os_ == s else 'FAIL'}")
    print(f"\ntotals: {OS}/{S} blocks, {OE}/{E} elements   ({time.time()-t0:.0f}s)")
    print("GREEDY-FROM-TOP HOLDS AT EVERY INTERMEDIATE STAGE" if OS == S
          else "*** fails ***")
    sys.exit(0 if OS == S else 1)
