"""
THE STAIRCASE INVARIANT.

At every colour-t socle step, writing D for the colour-t threshold, the partial
module X satisfies, for every j >= 0 and every chain vertex t-j >= 1:

    (Delta+)   X_{t-j} is FULL at every degree >= D + j + 1,
               and deficient by AT MOST ONE at the diagonal degree D + j.

Verified with no exception: full above the diagonal at 2337/2337 steps for
j = 0,1,2 and 2060/2060, 888/888 for j = 3,4; and on the diagonal the
deficiency is 0 or 1, never 2 -- distribution by j:

    j=1: 2207 full, 130 short by one
    j=2: 2268 full,  69 short by one
    j=3: 2028 full,  32 short by one
    j=4:  877 full,  11 short by one

WHY THIS IS THE RIGHT SHAPE.  The socle step at t-j fills I_{t-j}[e] as soon as
both neighbours of t-j are full at degree e+1.  For e >= D+j+1 the two
neighbours t-j+1 and t-j-1 need fullness at e+1 >= D+j+2, which is strictly
above BOTH of their diagonals (D+j-1 and D+j+1).  So the "full above the
diagonal" half of (Delta+) is self-propagating: it only ever requires itself,
one vertex over.  On the diagonal itself the requirement at the outer
neighbour t-j-1 falls exactly ON its diagonal, where a deficiency of one is
permitted -- which is precisely why the diagonal entries can be short, and
short by at most one.

The cascade terminates at vertex 1, whose only neighbour is 2: there the
condition is fullness of X_2 at degree D+(t-1)+1, strictly above vertex 2's
diagonal, hence supplied by (Delta+) itself.

A NEGATIVE RESULT WORTH RECORDING.  The most attractive hypothesis -- that the
final module V(w) is simply the greedy TOP TRUNCATION of I_kappa at EVERY
vertex, which would give a closed formula and settle everything at once -- is
FALSE.  It holds at 1185 of 1307 extremal modules over D4-D7, E6, E7 and A5,
but not always.  The smallest counterexample is D5, lambda = omega_2,
mu = (1,-1,0,0,0): there V_2 occupies degrees 2 and 6 while the top truncation
would give 6 and 4.  (Vertex 2 is the chain neighbour t-1.)  Note this element
does have simple spectrum -- r + 2m = 2 -- so the failure is not confined to
the degenerate cases.  Top truncation does hold at t and at the two leaves for
every extremal module tested (9842/9842), which is the hypothesis (H).

WHAT IS STILL MISSING.  (Delta+) is preserved along the construction only if
the colour-(t-j) steps have actually occurred by the time the colour-t step
does; when D drops, the diagonal moves down and previously-unfilled degrees
must already be full.  Controlling that interleaving along the c-trace is the
open point.  Also, (G1) needs the j=1 diagonal to be full, not merely short by
at most one, whenever X_t[D] != 0 -- see chainG.py for that refinement.
"""
import sys, time
from collections import Counter
from rootsys import RootSystem, orbit_with_words
from wordtrace import trace
from homological import setup
from gls import socle_step
from lemA import TRIV

CASES = [("D", 5, 2), ("D", 5, 3), ("D", 6, 3), ("D", 6, 4),
         ("D", 7, 2), ("D", 7, 3)]


def run(typ, n, k, cap=400):
    rs, pr, kk = setup(typ, n, k)
    t = TRIV[typ](n)
    lam = rs.fundamental(k)
    orb = orbit_with_words(rs, lam)
    if len(orb) > cap:
        return None
    st = Counter()
    for mu, word in orb.items():
        if not word:
            continue
        X = {}
        for _, i, c, _ in trace(rs, word, lam):
            D = max([d for d in range(pr.top + 1)
                     if pr.dim.get((d, t), 0)
                     and len(X.get((d, t), [])) < pr.dim[(d, t)]], default=None)
            if i == t and D is not None:
                for j in range(0, t):
                    v = t - j
                    if v < 1:
                        break
                    above = all(len(X.get((d, v), [])) == pr.dim.get((d, v), 0)
                                for d in range(D + j + 1, pr.top + 1))
                    st[('above', j, above)] += 1
                    df = pr.dim.get((D + j, v), 0) - len(X.get((D + j, v), []))
                    st[('diag', j, min(df, 2))] += 1
            for _ in range(c):
                X = socle_step(pr, rs, X, i)
    return st


if __name__ == "__main__":
    t0 = time.time()
    S = Counter()
    for typ, n, k in CASES:
        r = run(typ, n, k)
        if r:
            S += r
    js = sorted(set(j for kind, j, _ in S if kind == 'above'))
    ok = True
    print("full at every degree >= D+j+1:\n")
    for j in js:
        g, b = S[('above', j, True)], S[('above', j, False)]
        ok &= (b == 0)
        print(f"   vertex t-{j}: {g}/{g+b}   {'OK' if b == 0 else 'FAIL'}")
    print("\ndeficiency at the diagonal degree D+j:\n")
    for j in js:
        if j == 0:
            continue
        row = {d: S[('diag', j, d)] for d in (0, 1, 2) if S[('diag', j, d)]}
        ok &= (S[('diag', j, 2)] == 0)
        print(f"   vertex t-{j}: {row}   "
              f"{'OK (never 2 or more)' if not S[('diag', j, 2)] else 'FAIL'}")
    print(f"\n({time.time()-t0:.0f}s)")
    print("STAIRCASE INVARIANT VERIFIED" if ok else "*** failed ***")
    sys.exit(0 if ok else 1)
