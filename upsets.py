"""
THE CONJECTURE REDUCED TO ONE STATEMENT AT ONE VERTEX.

Carry along the socle construction the hypothesis

        (H)   X is an up-set at t and at the two leaves.

Steps at colours other than t and the leaves do not change X there, so (H) has
to be checked only at colour-t steps and at leaf steps.

LEAF STEPS -- the mechanism is proved.  At a colour-j step with j a leaf, whose
only neighbour is t, c_j = v'_t - 2 v'_j >= 1 because the word is reduced.  The
leaf degrees are g_1+1 > ... > g_k+1 with g_i = t+k-2i, so an up-set at j
occupies the top v'_j of them and the largest unfilled one is

        D_j = t + k - 2 v'_j - 1.

The colour-t threshold is D = t+k-2a if v'_t = 2a and D = t+k-2a-2 if
v'_t = 2a+1.  In the first case c_j >= 1 forces v'_j <= a-1, in the second
v'_j <= a; either way D_j > D.  Since X is an up-set at t it is FULL above D,
hence full at D_j + 1, so the socle condition  a_{j->t} xi in X_t[D_j+1]  is
automatic and the top unfilled leaf slot lies in soc_j(I/X).  The same
inequality applies to every unfilled leaf slot at a degree >= D, so the socle
contains that entire top block.

COLOUR-t STEPS.  Preservation of the up-set at t is (G1)/(G2), and both are now
proved from (H) plus the reducedness bound c >= 1 and two inequalities about
the profile of I_kappa alone:

    (G1)  cum(I_{t-1}, deg >= D+1) <= 2a   - [k=t]      (52/52)
    (G2)  cum(I_{t-1}, deg >= D+1) <= 2a+2 - [k=t]      (80/80)

together with, at the leaves, min(v'_{n-1},v'_n) >= a in the (G1) case, and in
the (G2) case exactly one leaf deficient by exactly one dimension when r' = 0
and none when r' = 2 -- all from
    v'_{n-1}+v'_n = v'_t - mu'_{n-1},  |v'_{n-1}-v'_n| = |mu'_n| <= 1,
and mu'_{n-1} <= 0, which itself follows from c = mu'_t - mu'_{n-1} >= 1 and
mu'_t <= 1.

WHAT IS LEFT.  Both (G1) and (G2) also need X to be an up-set at the CHAIN
neighbour t-1, which is not part of (H).  So the whole conjecture now rests on:

        (*)  at every colour-t socle step, X is an up-set at t-1.

Verified 2179/2179.  It is not part of a uniform statement: X is an up-set at
every neighbour of the colour being applied only 18982 times out of 19076, the
failures all occurring at colour-1 steps far from t; and globally X is an
up-set at t-1 at 8959 of 8960 stages, the single exception occurring at a
colour-1 step, where it is not needed.

This script verifies (H)'s preservation, the leaf mechanism, and (*).
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


def upset(pr, X, j):
    degs = [d for d in range(pr.top + 1) if pr.dim.get((d, j), 0)]
    prof = [(d, len(X.get((d, j), [])), pr.dim[(d, j)]) for d in degs]
    part = [d for d, a, b in prof if 0 < a < b]
    zero = [d for d, a, b in prof if a == 0]
    full = [d for d, a, b in prof if a == b]
    return len(part) <= 1 and not (zero and full and max(zero) > min(full))


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
            if i in (n - 1, n):
                jd = [d for d in range(pr.top + 1) if pr.dim.get((d, i), 0)]
                unf = [d for d in jd if len(X.get((d, i), [])) < pr.dim[(d, i)]]
                if unf and D is not None:
                    st['Dj_tot'] += 1
                    st['Dj_ok'] += (max(unf) > D)
            if i == t:
                st['star_tot'] += 1
                st['star_ok'] += upset(pr, X, t - 1)
            for _ in range(c):
                X = socle_step(pr, rs, X, i)
            st['H_tot'] += 1
            st['H_ok'] += (upset(pr, X, t) and upset(pr, X, n - 1)
                           and upset(pr, X, n))
    return st


if __name__ == "__main__":
    t0 = time.time()
    S = Counter()
    for typ, n, k in CASES:
        r = run(typ, n, k)
        if r:
            S += r
    def line(lbl, a, b):
        print(f"  {lbl:56} {S[a]}/{S[b]}   {'OK' if S[a] == S[b] else 'FAIL'}")
    line("(H) up-set at t and both leaves, after every step", 'H_ok', 'H_tot')
    line("leaf mechanism: D_j > D at every leaf step", 'Dj_ok', 'Dj_tot')
    line("(*) up-set at t-1 at every colour-t step", 'star_ok', 'star_tot')
    ok = all(S[a] == S[b] for a, b in [('H_ok', 'H_tot'), ('Dj_ok', 'Dj_tot'),
                                       ('star_ok', 'star_tot')])
    print(f"\n({time.time()-t0:.0f}s)")
    print("VERIFIED -- the conjecture rests on (*) alone" if ok
          else "*** something failed ***")
    sys.exit(0 if ok else 1)
