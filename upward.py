"""
LEMMA U (upward closure at the trivalent node, type D).

Let t be the trivalent node, with neighbours j_1,j_2,j_3.  For a homogeneous
x in I_t[d] the submodule Pi x contains the three elements

        y_r = a_{j_r -> t} a_{t -> j_r} x   in  I_t[d+2],

and the preprojective relation at t says  sum_r eps_r y_r = 0.  Hence
dim span(y_1,y_2,y_3) <= 2.  In type D, Lemma P gives dim I_t[d+2] <= 2, so
there is room for the span to be everything -- and this script checks that it
IS everything, for every nonzero x.

Consequence: if V is a graded submodule and V_t[d] != 0 then
V_t[d+2] = I_t[d+2]; inductively V_t is UPWARD CLOSED, and a submodule with
prescribed colour-t total v_t therefore has the greedy top-truncated profile.
That is Proposition T, for all k.

In type E, dim I_t[d] reaches 6, so the span of three elements subject to one
relation cannot fill it, and the argument -- correctly -- breaks down.
"""
import sys, time, itertools
from rootsys import RootSystem
from proj2 import Proj2
from qvar import rref, P
from homological import setup
from lemA import TRIV


def span_of_double_paths(pr, rs, t, d, x):
    """images a_{j->t} a_{t->j} x in I_t[d+2], for each neighbour j."""
    out = []
    for j in sorted(rs.adj[t]):
        if (d, t, j) not in pr.arrow:
            continue
        y = pr.act(d, t, x, j)
        if (d + 1, j, t) not in pr.arrow:
            continue
        out.append(pr.act(d + 1, j, y, t))
    return out


def run(typ, n, k):
    rs, pr, kk = setup(typ, n, k)
    t = TRIV[typ](n)
    tot = ok = 0
    worst = []
    for d in range(pr.top + 1):
        m = pr.dim.get((d, t), 0)
        m2 = pr.dim.get((d + 2, t), 0)
        if m == 0 or m2 == 0:
            continue
        for coeffs in itertools.product(range(P), repeat=m):
            if not any(coeffs):
                continue
            x = list(coeffs)
            ys = span_of_double_paths(pr, rs, t, d, x)
            R, _ = rref([y[:] for y in ys], m2)
            tot += 1
            if len(R) == m2:
                ok += 1
            elif len(worst) < 3:
                worst.append((d, m, m2, len(R)))
    return tot, ok, worst


CASES = [("D", 4, 1), ("D", 4, 2), ("D", 5, 2), ("D", 5, 3), ("D", 6, 2),
         ("D", 6, 3), ("D", 6, 4), ("D", 7, 3), ("D", 7, 5),
         ("E", 6, 3), ("E", 7, 4)]

if __name__ == "__main__":
    t0 = time.time()
    T = O = 0
    dfail = efail = 0
    for typ, n, k in CASES:
        tot, ok, worst = run(typ, n, k)
        if typ == "D":
            T += tot; O += ok
            if ok != tot:
                dfail += 1
        else:
            efail += (ok != tot)
        tag = "OK" if ok == tot else f"incomplete {worst}"
        print(f"  {typ}{n} om{k}: {ok}/{tot} vectors x with "
              f"span(double paths) = I_t[d+2]   {tag}")
    print(f"\ntype D totals: {O}/{T}   ({time.time()-t0:.0f}s)")
    print("LEMMA U HOLDS THROUGHOUT TYPE D" if O == T and dfail == 0
          else "*** fails in type D ***")
    print("(type E rows are expected to fail: dim I_t[d] there exceeds 2)")
    sys.exit(0 if O == T else 1)
