"""
LEMMA H.  Let X = soc_{(j_1,...,j_p)}(I) be an iterated socle with each
consecutive pair of colours distinct.  Then the head of X is concentrated at
the LAST colour j_p:

        Hom(X, S_l) = 0   for every  l != j_p.

Consequence, with the homological identity
   dim soc_i(I/X) = <alpha_i^vee, mu'> + dim top_i(X),
this gives  dim soc_i(I/X) = c  whenever i != j_p, which is exactly Lemma S:
the first of the c repeated steps adds all c dimensions at once, and each
further repeat adds  (-c) + c = 0.  So the c-weighted sequence is idempotent
on repeats and the total growth is c.

This script measures head(X) at every step of every element.
"""
import sys, time
from collections import Counter
from rootsys import RootSystem, orbit_with_words
from wordtrace import trace
from proj2 import Proj2
from gls import socle_step
from homological import setup, top_i


def run(typ, n, k, cap=400):
    rs, pr, kk = setup(typ, n, k)
    lam = rs.fundamental(k)
    orb = orbit_with_words(rs, lam)
    if len(orb) > cap:
        return None
    tested = ok = 0
    firststep_ok = firststep = 0
    ctr = Counter()
    for mu, word in orb.items():
        if not word:
            continue
        last = None
        X = {}
        for _, i, c, _ in trace(rs, word, lam):
            # --- Lemma H on the current X (built from colours ending in `last`)
            if last is not None:
                tested += 1
                heads = {l: top_i(pr, rs, X, l) for l in rs.I}
                nz = {l for l, h in heads.items() if h}
                ok += (nz == {last})
                if nz != {last}:
                    ctr[(typ, n, k, last, tuple(sorted(nz)))] += 1
            # --- does the FIRST repeat already add all c?
            before = sum(len(X.get((d, i), [])) for d in range(pr.top + 1))
            X = socle_step(pr, rs, X, i)
            after = sum(len(X.get((d, i), [])) for d in range(pr.top + 1))
            firststep += 1
            firststep_ok += (after - before == c)
            for _ in range(c - 1):
                X = socle_step(pr, rs, X, i)
            last = i
    return tested, ok, firststep, firststep_ok, ctr


CASES = [("A", 4, 2), ("A", 5, 3), ("D", 4, 2), ("D", 5, 2), ("D", 5, 3),
         ("D", 6, 3), ("E", 6, 2), ("E", 6, 3), ("E", 6, 5), ("E", 7, 1)]

if __name__ == "__main__":
    t0 = time.time()
    T = O = F = FO = 0
    C = Counter()
    for typ, n, k in CASES:
        r = run(typ, n, k)
        if r is None:
            continue
        t, o, f, fo, c = r
        T += t; O += o; F += f; FO += fo; C += c
        print(f"  {typ}{n} om{k}: head concentrated at last colour {o}/{t}; "
              f"first repeat adds all c {fo}/{f}")
    print(f"\ntotals: Lemma H {O}/{T},  one-step-suffices {FO}/{F}  "
          f"({time.time()-t0:.0f}s)")
    if C:
        print("deviations (type,n,k,last colour,head support):")
        for kk_, vv in C.most_common(10):
            print("   ", kk_, vv)
    sys.exit(0 if O == T and FO == F else 1)
