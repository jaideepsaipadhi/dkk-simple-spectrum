"""
LEMMA S (the dimension count).

    Let lambda = omega_kappa, w in W, mu = w.lambda, and let sigma be the
    c-weighted sequence of a reduced word for w (each c-trace step with letter
    i and value c contributes i repeated c times, in trace order).  Put
    I_kappa = Pi e_{nu^{-1}(kappa)} and X = soc_sigma(I_kappa).  Then

            dim X  =  v  =  lambda - mu   (in root coordinates),

    and moreover each individual socle step of colour i and value c raises
    dim X_i by exactly c.

Why this settles the construction:  by Savage-Tingley (Prop 4.9) a w-extremal
dimension vector v is carried by a UNIQUE submodule of I_kappa.  So any
submodule with dim = v IS the fixed-point module V.  Lemma S exhibits one.

This script tests both halves across A, D, E.
"""
import sys, time
from rootsys import RootSystem, orbit_with_words
from wordtrace import trace, v_from_trace
from proj2 import Proj2
from gls import socle_step


def socv(rs, j):
    pr = Proj2(rs, j)
    return max((i for (d, i), m in pr.dim.items() if d == pr.top and m), default=None)


def setup(typ, n, k, _c={}):
    key = (typ, n, k)
    if key in _c:
        return _c[key]
    rs = RootSystem(typ, n)
    nu = {j: socv(rs, j) for j in rs.I}
    inv = {b: a for a, b in nu.items()}
    _c[key] = (rs, Proj2(rs, inv[k]))
    return _c[key]


def dims(rs, X):
    return tuple(sum(len(X.get((d, i), [])) for d in range(60)) for i in rs.I)


def run(typ, n, k, cap=400):
    rs, pr = setup(typ, n, k)
    lam = rs.fundamental(k)
    orb = orbit_with_words(rs, lam)
    if len(orb) > cap:
        return None
    tot = okv = okstep = steps = 0
    for mu, word in orb.items():
        if not word:
            continue
        tr = trace(rs, word, lam)
        v = tuple(int(x) for x in v_from_trace(rs, tr))
        X = {}
        prev = dims(rs, X)
        good = True
        for _, i, c, _ in tr:
            for _ in range(c):
                X = socle_step(pr, rs, X, i)
            cur = dims(rs, X)
            steps += 1
            grew = cur[i - 1] - prev[i - 1]
            other = [j for j in rs.I if j != i and cur[j - 1] != prev[j - 1]]
            if grew == c and not other:
                okstep += 1
            else:
                good = False
            prev = cur
        tot += 1
        okv += (dims(rs, X) == v)
    return tot, okv, steps, okstep


CASES = [("A", 4, 2), ("A", 5, 3),
         ("D", 4, 1), ("D", 4, 2), ("D", 4, 3),
         ("D", 5, 2), ("D", 5, 3), ("D", 5, 4),
         ("D", 6, 2), ("D", 6, 3),
         ("E", 6, 1), ("E", 6, 2), ("E", 6, 3), ("E", 6, 5), ("E", 6, 6),
         ("E", 7, 1), ("E", 7, 7)]

if __name__ == "__main__":
    t0 = time.time()
    T = OKV = S = OKS = 0
    allok = True
    for typ, n, k in CASES:
        r = run(typ, n, k)
        if r is None:
            print(f"  {typ}{n} om{k}: skipped (orbit too large)")
            continue
        tot, okv, steps, okstep = r
        T += tot; OKV += okv; S += steps; OKS += okstep
        good = (okv == tot and okstep == steps)
        allok &= good
        print(f"  {typ}{n} om{k}: dim X = v in {okv}/{tot} elements; "
              f"per-step growth exact in {okstep}/{steps} steps  "
              f"{'OK' if good else 'FAIL'}")
    print(f"\ntotals: {OKV}/{T} elements, {OKS}/{S} socle steps  ({time.time()-t0:.0f}s)")
    print("LEMMA S HOLDS ON EVERY CASE TESTED" if allok else "*** LEMMA S FAILS ***")
    sys.exit(0 if allok else 1)
