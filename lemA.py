"""
LEMMA A, tested on the fixed-point modules V themselves (not on the ambient
injective).  Two candidate forms:

  (A1)  max_d dim V_i[d]  <=  max_d dim V_t[d]   for every vertex i
        (t = the trivalent node).  This implies "t simple => all simple".

  (A2)  the weaker statement itself:  V_t multiplicity-free  =>  V
        multiplicity-free.

V is computed by Theorem A (the c-weighted iterated socle), which is proved.
"""
import sys, time
from collections import Counter
from rootsys import RootSystem, orbit_with_words
from wordtrace import trace
from proj2 import Proj2
from gls import socle_step
from homological import setup

TRIV = {"D": lambda n: n - 2, "E": lambda n: 4, "A": lambda n: None}


def build(rs, pr, word, lam):
    X = {}
    for _, i, c, _ in trace(rs, word, lam):
        for _ in range(c):
            X = socle_step(pr, rs, X, i)
    return X


def run(typ, n, k, cap=400):
    rs, pr, kk = setup(typ, n, k)
    t = TRIV[typ](n)
    lam = rs.fundamental(k)
    orb = orbit_with_words(rs, lam)
    if len(orb) > cap:
        return None
    n1 = n2 = tot = 0
    bad = Counter()
    for mu, word in orb.items():
        if not word:
            continue
        X = build(rs, pr, word, lam)
        mx = {i: max((len(X.get((d, i), [])) for d in range(pr.top + 1)),
                     default=0) for i in rs.I}
        tot += 1
        n1 += all(mx[i] <= mx[t] for i in rs.I)
        n2 += (mx[t] <= 1) == (max(mx.values()) <= 1)
        if not all(mx[i] <= mx[t] for i in rs.I):
            bad[(typ, n, k, tuple(mx[i] for i in rs.I))] += 1
    return tot, n1, n2, bad


CASES = [("D", 4, 2), ("D", 5, 2), ("D", 5, 3), ("D", 6, 2), ("D", 6, 3),
         ("E", 6, 2), ("E", 6, 3), ("E", 6, 5), ("E", 7, 1), ("E", 7, 7),
         ("E", 7, 2), ("E", 7, 6)]

if __name__ == "__main__":
    t0 = time.time()
    T = N1 = N2 = 0
    B = Counter()
    for typ, n, k in CASES:
        r = run(typ, n, k)
        if r is None:
            print(f"  {typ}{n} om{k}: skipped")
            continue
        tot, a1, a2, b = r
        T += tot; N1 += a1; N2 += a2; B += b
        print(f"  {typ}{n} om{k}: (A1) max_i <= max_t {a1}/{tot};  "
              f"(A2) t simple => all simple {a2}/{tot}")
    print(f"\ntotals: (A1) {N1}/{T}   (A2) {N2}/{T}   ({time.time()-t0:.0f}s)")
    if B:
        print("(A1) failures — profile of max multiplicities by vertex:")
        for kk_, vv in B.most_common(8):
            print("   ", kk_, vv)
    sys.exit(0 if N2 == T else 1)
