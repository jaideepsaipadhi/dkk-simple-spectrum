"""
THE CHAIN NEIGHBOUR AT DEGREE D+1 -- what is proved, and what is left.

The conjecture reduces (see upsets.py) to a statement about ONE graded piece:

    (P1)  at a colour-t socle step, if X_t[D] != 0 then X_{t-1}[D+1] is full,
          and if X_t[D] = 0 then X_{t-1}[D+1] is deficient by at most one.

(P1) is strictly weaker than the up-set property at t-1: X_{t-1} is full in all
degrees > D at only 2207 of 2337 colour-t steps, whereas (P1) holds at all of
them.  So degree D+1 is exactly the right place to look.

PROVED.  Suppose X_t[D] != 0, so it is a line spanned by a SPECIAL vector x --
annihilated by an arrow from t to one of the leaves, say a_{t->n}.  Then

        a_{t->(t-1)} x  !=  0.

For if it vanished, x would be killed by two of the three arrows out of t; it
cannot be killed by the third as well, since I_kappa has simple socle S_kappa
and kappa <= n-2 is not a leaf, so y = a_{t->(n-1)} x is nonzero; but then the
preprojective relation at t reads

    eps_1 a_{(t-1)->t} a_{t->(t-1)} x + eps_2 a_{(n-1)->t} y
        + eps_3 a_{n->t} a_{t->n} x = 0,

whose first and third terms vanish, giving a_{(n-1)->t} y = 0.  Then y is a
nonzero element of the leaf n-1 killed by its only arrow, hence lies in
soc(I_kappa) = S_kappa, forcing kappa to be that leaf -- excluded.

So a_{t->(t-1)}(X_t[D]) is a nonzero line inside X_{t-1}[D+1], and (P1) holds
whenever dim I_{t-1}[D+1] = 1.

WHAT IS LEFT, exactly.  Splitting the colour-t steps by branch and by
r' = |supp mu' cap {n-1,n}| gives a completely clean table (2337 steps):

    branch            n     full   deficient by 1   by more
    X_t[D] != 0      547     547         0            0
    X_t[D] = 0, r'=0 855     855         0            0
    X_t[D] = 0, r'=1 935     805       130            0

The first row is what (G1) needs and the second and third are what (G2) needs:
deficiency never exceeds one, and it occurs only in the r'=1 branch, which is
exactly the branch where both leaves are full, so at most one neighbour is ever
deficient.  The three statements still to prove are therefore

    (A)  row 1 with dim I_{t-1}[D+1] = 2   (125 of the 547; the other 422 are
         proved by the nonvanishing argument above);
    (B)  row 2, always full;
    (C)  row 3, deficiency at most one -- automatic when dim I_{t-1}[D+1] = 1,
         so only the 2-dimensional case bites.

All three are about getting dimensions into the single space X_{t-1}[D+1].  The
natural mechanism is the socle step at t-1: I_{t-1}[D+1] enters soc_{t-1}(I/X)
as soon as X_t[D+2] and X_{t-2}[D+2] are full, and X_t[D+2] IS full because D
is the colour-t threshold.  So everything turns on X_{t-2}[D+2], one vertex
further out -- true at 2268 of 2337 steps -- and the cascade continues down the
chain, terminating when D+j passes the top degree of vertex t-j.  Making that
precise needs control of how colour-(t-j) steps interleave with colour-t steps,
which is where this argument currently stops.

This script verifies the nonvanishing, the split by dim I_{t-1}[D+1], and (P1).
"""
import sys, time
from collections import Counter
from rootsys import RootSystem, orbit_with_words
from wordtrace import trace
from homological import setup
from gls import socle_step
from lemA import TRIV
from qvar import rref, P

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
                m = pr.dim.get((D + 1, t - 1), 0)
                have = len(X.get((D + 1, t - 1), []))
                if X.get((D, t)):
                    x = X[(D, t)][0]
                    y = pr.act(D, t, x, t - 1) if (D, t, t - 1) in pr.arrow else []
                    st['nz_tot'] += 1
                    st['nz_ok'] += any(v % P for v in y)
                    st['dim%d' % m] += 1
                    st['P1_tot'] += 1
                    st['P1_ok'] += (have == m)
                else:
                    st['P1_tot'] += 1
                    st['P1_ok'] += (m - have <= 1)
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
    def line(lbl, a, b):
        print(f"  {lbl:54} {S[a]}/{S[b]}   {'OK' if S[a] == S[b] else 'FAIL'}")
    line("a_(t->t-1) nonzero on the special line X_t[D]", 'nz_ok', 'nz_tot')
    print(f"  of those, dim I_(t-1)[D+1] = 1 (PROVED case)         {S['dim1']}")
    print(f"                            = 2 (open case)            {S['dim2']}")
    line("(P1) at every colour-t step", 'P1_ok', 'P1_tot')
    ok = S['nz_ok'] == S['nz_tot'] and S['P1_ok'] == S['P1_tot']
    print(f"\n({time.time()-t0:.0f}s)")
    print("VERIFIED" if ok else "*** something failed ***")
    sys.exit(0 if ok else 1)
