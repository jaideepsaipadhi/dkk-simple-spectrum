"""
THE MISSING STATEMENT IN EPSILON-COORDINATES.

Everything below is type D_n with lambda = omega_k, k <= n-2, t = n-2.

In Bourbaki coordinates alpha_t = e_{n-2} - e_{n-1}, so for the intermediate
weight mu' = lambda - dim X of the c-trace,

        c = <alpha_t^vee, mu'> = mu'_{n-2} - mu'_{n-1},

and since mu' has entries in {0,+-1} and c >= 1 for a reduced word,

        c in {1,2},   and   c = 2  <=>  (mu'_{n-2}, mu'_{n-1}) = (1,-1).

Write r' = |supp mu' cap {n-1,n}| and m' = #{minus signs among the first n-2
coordinates}, so v'_t = r' + 2m' by Lemma C.

TWO STEPS THAT ARE PROVED.

  1. A critical step (one whose socle reaches degree D-2) has X_t[D] != 0.
     Indeed if X_t[D] = 0 then the double paths of the socle element at D-2 all
     vanish, and the propagation argument forces kappa to be a leaf or kappa=1,
     both excluded.

  2. Granting top-truncation for X, the greedy top truncation of the profile
     (1, 2^{k-1}, 1) to total v'_t has X_t[D] != 0 exactly when v'_t is EVEN and
     at least 2.  Hence a critical step has r' + 2m' even, i.e. r' in {0,2}.
     And r' = 0 forces mu'_{n-1} = 0, whence c = mu'_{n-2} <= 1, so c = 1.

So a critical step must have r' = 2.  Among steps with r' = 2 there are exactly
two signatures, (mu'_{n-2}, mu'_{n-1}) = (1,-1) with c = 2 and (0,-1) with
c = 1, and the whole remaining gap is:

        WHEN r' = 2 AND mu'_{n-2} = 0 (so c = 1), THE SOCLE STEP LANDS AT D,
        NOT AT D-2.

This script tabulates every colour-t step by (r', c, mu'_{n-2}, mu'_{n-1}) and
splits by whether the step is critical.  Observed: the critical steps form a
SINGLE class, (2, 2, 1, -1).
"""
import sys, time
from collections import Counter
from rootsys import RootSystem, orbit_with_words
from wordtrace import trace
from homological import setup
from gls import socle_step
from lemA import TRIV
from epsilon import eps_of_mu

CASES = [("D", 5, 2), ("D", 5, 3), ("D", 6, 2), ("D", 6, 3), ("D", 6, 4),
         ("D", 7, 2), ("D", 7, 3)]


def run(typ, n, k, cap=400):
    rs, pr, kk = setup(typ, n, k)
    t = TRIV[typ](n)
    lam = rs.fundamental(k)
    orb = orbit_with_words(rs, lam)
    if len(orb) > cap:
        return None
    crit, noncrit = Counter(), Counter()
    for mu, word in orb.items():
        if not word:
            continue
        X = {}
        for _, i, c, _ in trace(rs, word, lam):
            if i != t:
                for _ in range(c):
                    X = socle_step(pr, rs, X, i)
                continue
            vp = [sum(len(X.get((d, j), [])) for d in range(pr.top + 1))
                  for j in rs.I]
            e = eps_of_mu(n, k, vp)
            rp = sum(1 for u in (n - 2, n - 1) if e[u] != 0)
            D = max([d for d in range(pr.top + 1)
                     if pr.dim.get((d, t), 0)
                     and len(X.get((d, t), [])) < pr.dim[(d, t)]], default=None)
            Y = socle_step(pr, rs, X, t)
            iscrit = (D is not None
                      and len(Y.get((D - 2, t), [])) > len(X.get((D - 2, t), [])))
            key = (rp, c, e[n - 3], e[n - 2])
            (crit if iscrit else noncrit)[key] += 1
            X = Y
            for _ in range(c - 1):
                X = socle_step(pr, rs, X, t)
    return crit, noncrit


if __name__ == "__main__":
    t0 = time.time()
    C, N = Counter(), Counter()
    for typ, n, k in CASES:
        r = run(typ, n, k)
        if r is None:
            continue
        c, nc = r
        C += c; N += nc
    print("colour-t steps keyed by (r', c, mu'_{n-2}, mu'_{n-1})\n")
    print("  CRITICAL (socle reaches D-2):")
    for key, v in sorted(C.items()):
        print(f"      {key}   x{v}")
    print("  non-critical:")
    for key, v in sorted(N.items()):
        print(f"      {key}   x{v}")
    ok = set(C) == {(2, 2, 1, -1)}
    print(f"\ncritical steps form the single class (2,2,1,-1): "
          f"{'YES' if ok else 'NO -- ' + str(set(C))}")
    # the two proved implications
    bad_par = [key for key in C if (key[0] % 2) != 0]
    bad_r0 = [key for key in C if key[0] == 0]
    print(f"  every critical step has r' even   : {not bad_par}")
    print(f"  no critical step has r' = 0       : {not bad_r0}")
    print(f"  remaining gap: steps with r'=2, mu'_{{n-2}}=0 (c=1) are all "
          f"non-critical: "
          f"{all(key != (2, 1, 0, -1) for key in C)}")
    print(f"\n({time.time()-t0:.0f}s)")
    sys.exit(0 if ok else 1)
