"""
Proposition "top truncation": for V = soc_sigma(I_kappa) (a SUBMODULE of
I_kappa of dimension vector v = lambda - w.lambda), the colour-t multiplicity
sequence of V is the greedy TOP truncation of the colour-t sequence of
I_kappa to total v_t -- the same degrees, filled from the highest downward.
"""
import sys, time
from rootsys import RootSystem, orbit_with_words
from wordtrace import trace
from proj2 import Proj2
from homological import setup
from lemA import build, TRIV


def top_truncate(seq_by_deg, total):
    """seq_by_deg: list of (degree, multiplicity) in increasing degree."""
    out = {}
    left = total
    for d, m in reversed(seq_by_deg):
        take = min(m, left)
        if take:
            out[d] = take
        left -= take
        if left == 0:
            break
    return out, left


def run(typ, n, k, cap=400):
    rs, pr, kk = setup(typ, n, k)
    t = TRIV[typ](n)
    lam = rs.fundamental(k)
    orb = orbit_with_words(rs, lam)
    if len(orb) > cap:
        return None
    prof = [(d, pr.dim[(d, t)]) for d in range(pr.top + 1)
            if pr.dim.get((d, t), 0)]
    tot = ok = 0
    for mu, word in orb.items():
        if not word:
            continue
        X = build(rs, pr, word, lam)
        got = {d: len(X[(d, t)]) for d in range(pr.top + 1) if X.get((d, t))}
        vt = sum(got.values())
        want, left = top_truncate(prof, vt)
        tot += 1
        ok += (got == want and left == 0)
    return tot, ok, [m for _, m in prof]


if __name__ == "__main__":
    t0 = time.time()
    T = O = 0
    for typ, n in [("D", 4), ("D", 5), ("D", 6), ("D", 7)]:
        for k in range(1, n - 1):
            r = run(typ, n, k)
            if r is None:
                print(f"  {typ}{n} om{k}: skipped (orbit too large)")
                continue
            tot, ok, prof = r
            T += tot; O += ok
            print(f"  {typ}{n} om{k}: {ok}/{tot}   (I_k colour-t profile {prof})"
                  f"   {'OK' if ok == tot else 'FAIL'}")
    print(f"\ntotal {O}/{T}   ({time.time()-t0:.0f}s)")
    print("TOP TRUNCATION HOLDS" if O == T else "*** FAILS ***")
    sys.exit(0 if O == T else 1)
