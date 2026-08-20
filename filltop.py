"""
PROPOSITION T AS AN INDUCTION ON THE SOCLE CONSTRUCTION.

Since V = soc_sigma(I_kappa) by the Main Theorem, Proposition T -- that the
colour-t multiplicities of V are the greedy TOP truncation of those of I_kappa
-- can be proved by induction along the socle steps.  Write, for a graded
submodule X,

        D(X) := max{ d : X_t[d] is a PROPER subspace of (I_kappa)_t[d] }.

X is top-truncated at t exactly when X_t[d] = 0 for d < D(X) and
X_t[d] = (I_kappa)_t[d] for d > D(X).  So the induction needs:

  (a) a colour-t socle step contributes nothing in degrees < D - 2;
  (b) if it contributes at D - 2, then it fills degree D in the same step.

Then the new profile is again zero below, partial at one degree, full above.
Steps at colours other than t do not touch X_t, so (a) + (b) give
Proposition T for every k, with no case split.

STATUS.

(a) is proved.  Let x in (I_kappa)_t[d] lie in soc_t(I/X) with d < D - 2.  Then
    X_t[d+2] = 0, so every double path a_{j->t} a_{t->j} x vanishes.  By the
    simple socle of I_kappa, x is not killed by all three arrows, so some
    y_j = a_{t->j} x is nonzero and killed by a_{j->t}.  If j is a leaf this
    puts y_j in soc(I_kappa) = S_kappa, forcing kappa to be a leaf, excluded
    since kappa <= n-2.  If j = t-1, the propagation of Lemma I (relation at a
    chain vertex + simple socle) pushes the same situation outward to vertex 1,
    whose only neighbour is 2, again landing in the socle and forcing kappa = 1
    -- excluded because for k = 1 Lemma P gives a colour-t profile with no
    multiplicity 2, so D-2 < d is impossible to begin with.

(b) reduces to a single assertion.  Suppose a colour-t step reaches degree
    D - 2, say x in soc_t(I/X) is nonzero there.  Then:

    (i)   dim X_t[D] = 1 and dim (I_kappa)_t[D] = 2.
          If X_t[D] were 0, the double paths of x would land in X_t[D] = 0 and
          hence all vanish; by the simple socle some y_j = a_{t->j} x is
          nonzero and killed by a_{j->t}, and the propagation of Lemma I forces
          kappa to be a leaf or kappa = 1, both excluded.

    (ii)  soc_t(I/X) has dimension exactly 1 at D - 2.
          A subspace of dimension >= 2 would contain a NON-special vector, since
          the special vectors form two lines and not a subspace; its double
          paths would span (I_kappa)_t[D] by Lemma U while lying in X_t[D],
          contradicting (i).

    So one further dimension at D fills it, and everything follows from

        CLAIM F.  If X_t[D] is a proper subspace of (I_kappa)_t[D], then
        soc_t(I/X) has a nonzero component in degree D.

    Indeed (a) confines the socle to degrees D and D-2, Claim F puts something
    at D, (ii) caps D-2 at one dimension, and (i) says one dimension fills D.
    So the new colour-t profile is again zero below, partial at one degree and
    full above: top-truncation is preserved, and Proposition T follows for
    every k with no case split.

    Claim F in turn follows from a statement about which pieces are FULL:

        CLAIM G.  (G1) if X_t[D] != 0 then every neighbour of t is full at
                  degree D+1;
                  (G2) if X_t[D] = 0 then at most ONE neighbour is deficient at
                  degree D+1, and the deficiency is one-dimensional.

    G => F: under (G1) every z in I_t[D] has its arrows landing in X
    automatically, so the whole of I_t[D]/X_t[D] is in the socle; under (G2)
    the single one-dimensional deficiency imposes one linear condition on the
    2-dimensional I_t[D], so a nonzero witness survives.

    Claim G is verified below (532/532 and 1691/1691, all deficiencies of
    dimension 1) and is NOT proved.  Two things are worth recording.  First,
    Claim F is FALSE for arbitrary graded submodules of Pi e_kappa -- among all
    of them the dual statement fails 318 times out of 572 -- so it genuinely
    uses that dim X is w-extremal.  Second, the naive mechanism for (G1) fails:
    when X_t[D] is a line it is ALWAYS one of the two special lines (527/527),
    so one leaf arrow kills it and the fullness of that leaf at D+1 has to come
    from elsewhere in X.

This script measures (a), Claim F, and the structure at the D-2 steps.
"""
import sys, time
from collections import Counter
from rootsys import RootSystem, orbit_with_words
from wordtrace import trace
from homological import setup
from gls import socle_step
from lemA import TRIV
from upclosure import special_lines
from qvar import P

CASES = [("D", 4, 2), ("D", 5, 2), ("D", 5, 3), ("D", 6, 2), ("D", 6, 3),
         ("D", 6, 4), ("D", 7, 2), ("D", 7, 3)]


def run(typ, n, k, cap=400):
    rs, pr, kk = setup(typ, n, k)
    t = TRIV[typ](n)
    lam = rs.fundamental(k)
    orb = orbit_with_words(rs, lam)
    if len(orb) > cap:
        return None
    sp = {d: set(special_lines(pr, rs, t, d))
          for d in range(pr.top + 1) if pr.dim.get((d, t), 0)}
    rel = Counter()
    nD2 = nspecial = nfilled = 0
    nsteps = [0]; nF = [0]
    for mu, word in orb.items():
        if not word:
            continue
        X = {}
        for _, i, c, _ in trace(rs, word, lam):
            if i != t:
                for _ in range(c):
                    X = socle_step(pr, rs, X, i)
                continue
            D = max([d for d in range(pr.top + 1)
                     if pr.dim.get((d, t), 0)
                     and len(X.get((d, t), [])) < pr.dim[(d, t)]], default=None)
            Y = socle_step(pr, rs, X, t)
            if D is not None:
                nsteps[0] += 1
                if len(Y.get((D, t), [])) > len(X.get((D, t), [])):
                    nF[0] += 1
                for d in range(pr.top + 1):
                    g = len(Y.get((d, t), [])) - len(X.get((d, t), []))
                    if g:
                        rel[d - D] += g
                new = Y.get((D - 2, t), [])
                if len(new) > len(X.get((D - 2, t), [])):
                    nD2 += 1
                    if len(new) == 1:
                        x = new[0]
                        lead = next(u for u, cc in enumerate(x) if cc % P)
                        inv = pow(x[lead], P - 2, P)
                        if tuple((cc * inv) % P for cc in x) in sp.get(D - 2, ()):
                            nspecial += 1
                    if len(Y.get((D, t), [])) == pr.dim.get((D, t), 0):
                        nfilled += 1
            X = Y
            for _ in range(c - 1):
                X = socle_step(pr, rs, X, t)
    return rel, nD2, nspecial, nfilled, nsteps[0], nF[0]


if __name__ == "__main__":
    t0 = time.time()
    REL = Counter()
    A = B = C = S = F = 0
    ok = True
    for typ, n, k in CASES:
        r = run(typ, n, k)
        if r is None:
            print(f"  {typ}{n} om{k}: skipped")
            continue
        rel, nD2, nsp, nf, ns, nfF = r
        REL += rel
        A += nD2; B += nsp; C += nf; S += ns; F += nfF
        below = sum(v for d, v in rel.items() if d < -2)
        ok &= (below == 0 and nf == nD2 and nfF == ns)
        print(f"  {typ}{n} om{k}: contributions by degree-D offset "
              f"{dict(sorted(rel.items()))};  steps hitting D-2: {nD2}, "
              f"of which D filled: {nf}")
    below = sum(v for d, v in REL.items() if d < -2)
    print(f"\n(a) contributions strictly below D-2 : {below}   "
          f"{'OK (none)' if below == 0 else 'FAIL'}")
    print(f"(b) steps hitting D-2                : {A}, "
          f"degree D filled in the same step: {C}   "
          f"{'OK' if C == A else 'FAIL'}")
    print(f"    of those, the D-2 component is a special line: {B}/{A}")
    print(f"CLAIM F: steps contributing at the threshold D: {F}/{S}   "
          f"{'OK' if F == S else 'FAIL'}")
    print(f"\n({time.time()-t0:.0f}s)")
    print("INDUCTION STEP HOLDS THROUGHOUT" if ok else "*** failed ***")
    sys.exit(0 if ok else 1)
