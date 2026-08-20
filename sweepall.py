"""
Aggregate agreement between the criterion (Corollary: no repeated
(colour, degree) in soc_sigma(I_kappa)) and the INDEPENDENT solver dfs3.solve,
which searches graded submodules of the projective directly.

Prints one line per (type, rank, fundamental weight) and a grand total.
Set DKKP to test in a different characteristic, e.g.  DKKP=11 python sweepall.py
"""
import os, sys, time, json
from rootsys import RootSystem, orbit_with_words
from wordtrace import trace, v_from_trace
from dfs3 import solve
from homological import setup
from lemA import build

CASES = ([("A", n, k) for n in range(3, 7) for k in range(1, n + 1)] +
         [("D", n, k) for n in range(4, 8) for k in range(1, n + 1)] +
         [("E", 6, k) for k in range(1, 7)] +
         [("E", 7, k) for k in (1, 7)])

ORBIT_CAP = int(os.environ.get("DKKCAP", 400))


def main():
    t0 = time.time()
    T = A = 0
    per = {}
    for typ, n, k in CASES:
        rs0 = RootSystem(typ, n)
        orb = orbit_with_words(rs0, rs0.fundamental(k))
        if len(orb) > ORBIT_CAP:
            print(f"  {typ}{n} om{k}: skipped, orbit {len(orb)} > {ORBIT_CAP}")
            continue
        rs, pr, kk = setup(typ, n, k)
        lam = rs.fundamental(k)
        tot = agree = 0
        for mu, word in orb.items():
            if not word:
                continue
            v = tuple(int(x) for x in v_from_trace(rs, trace(rs, word, lam)))
            res = solve(rs, k, v)
            if not res:
                continue
            g = list(res.values())[0]
            truth = all(c <= 1 for row in g.values() for c in row.values())
            X = build(rs, pr, word, lam)
            pred = all(len(x) <= 1 for x in X.values())
            tot += 1
            agree += (pred == truth)
        T += tot; A += agree
        per[typ] = per.get(typ, [0, 0])
        per[typ][0] += agree; per[typ][1] += tot
        print(f"  {typ}{n} om{k}: {agree}/{tot}   {'OK' if agree==tot else 'FAIL'}")
    print()
    for typ in sorted(per):
        print(f"  type {typ}: {per[typ][0]}/{per[typ][1]}")
    print(f"\nTOTAL {A}/{T}   (p = {os.environ.get('DKKP','7')}, "
          f"{time.time()-t0:.0f}s)")
    json.dump({"total": T, "agree": A, "per": per}, open("sweepall.json", "w"))
    sys.exit(0 if A == T else 1)


if __name__ == "__main__":
    main()
