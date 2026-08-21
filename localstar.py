"""
IS (*o) LOCAL?  NO -- extremality cannot be dropped.

(*o) (starchar.py) says: an extremal module M <= I_kappa with
c_t = <alpha_t^vee, lambda - dim M> >= 1 is an up-set at the chain neighbour
t-1.  It is order-free, but it still refers to extremality.  The natural next
reduction would be to drop that too, keeping only the inductive hypothesis (H)
(up-set at t and at the two leaves), i.e.

    (L)   for EVERY graded submodule N <= I_kappa,
          (H) and c_t >= 1  ==>  N is an up-set at t-1.

(L) is FALSE.  It holds for every graded submodule of I_kappa with kappa <= 3,
which is what makes it tempting, but it fails at kappa = 4 -- where the chain
neighbour t-1 first carries multiplicity 2 in two different degrees.  The
smallest counterexample is D_6, kappa = 4, t = 4, the graded submodule with

    vertex 1  0/1 1/1
    vertex 2  0/1 1/2 1/1
    vertex 3  0/1 1/2 1/2 1/1        <- two partial degrees: not an up-set
    vertex 4  0/1 0/2 1/2 2/2 1/1    <- up-set at t
    vertex 5  0/1 1/1 1/1 1/1        <- up-set at the leaf
    vertex 6  0/1 0/1 1/1 1/1        <- up-set at the leaf

(entries are dim N / dim I_kappa by degree), which has c_t = 1.  There are 2
such submodules for D_6, kappa=4 and 2 for D_7, kappa=4; without (H) there are
10 and 10.

So (*o) genuinely needs its modules to be extremal, and the remaining gap does
not reduce to a statement about I_kappa alone.  This script records that.
"""
import sys, time
from collections import Counter
from rootsys import RootSystem
from homological import setup
from allsub import all_submodules
from upsets import upset

CASES = [("D", 4, 2), ("D", 5, 2), ("D", 5, 3), ("D", 6, 2),
         ("D", 6, 3), ("D", 6, 4), ("D", 7, 2), ("D", 7, 3), ("D", 7, 4)]


def run(typ, n, k):
    rs, pr, kk = setup(typ, n, k)
    t = n - 2
    _, subm = all_submodules(rs, kk)
    st = Counter()
    for N in subm:
        vv = {j: sum(len(N.get((d, j), [])) for d in range(pr.top + 1))
              for j in rs.I}
        ct = (1 if t == kk else 0) + sum(vv[j] for j in rs.adj[t]) - 2 * vv[t]
        H = upset(pr, N, t) and upset(pr, N, n - 1) and upset(pr, N, n)
        st['all'] += 1
        if ct >= 1:
            st['tot'] += 1
            st['ok'] += upset(pr, N, t - 1)
            if H:
                st['Htot'] += 1
                st['Hok'] += upset(pr, N, t - 1)
    return st


if __name__ == "__main__":
    t0 = time.time()
    S = Counter()
    for typ, n, k in CASES:
        r = run(typ, n, k)
        S += r
        flag = '' if r['Hok'] == r['Htot'] else '   <-- (L) FAILS'
        print(f"  {typ}{n} kappa={k}: (H) and c_t>=1 -> up-set at t-1 "
              f"{r['Hok']}/{r['Htot']}  (of {r['all']} graded submodules){flag}")
    print(f"\n  (L) over all graded submodules with (H): {S['Hok']}/{S['Htot']}")
    print(f"  without (H):                            {S['ok']}/{S['tot']}")
    ok = S['Hok'] < S['Htot']
    print(f"\n({time.time()-t0:.0f}s)")
    print("(L) IS FALSE, as expected -- extremality is needed" if ok
          else "*** (L) did not fail: the recorded counterexamples are gone ***")
    sys.exit(0 if ok else 1)
