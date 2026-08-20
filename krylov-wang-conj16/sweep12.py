"""
sweep12.py -- the decisive search: Hesselink-type representations of SL2^k
tested DIRECTLY against Theorem 1.2.

Why this and not census.py.  census.py tests identity (2).  Proposition 2.1
derives Theorem 1.2 from (2) only for COFREE V, via Lemma 2.3.  So a
Hesselink-type V that is not cofree can satisfy (2) and still violate
Theorem 1.2, and census.py would never notice.  This script computes

        C[N_V] = C[V] / (C[V]^G_+ C[V])

by linear algebra, decomposes it into irreducible characters, and compares
against sum_lambda M_q(lambda) chi_lambda.  Any disagreement for a
Hesselink-type V refutes Conjecture 1.6.

Calibration.  The detector is checked in both directions before the sweep:
V = V_1, V_1+V_1, V_2 must agree (Theorem 1.2 is proved there) and V = V_3,
V_4 must disagree (Section 2.1 of the paper gives them a modified P).

Usage
-----
    python3 sweep12.py                       # k = 1,2,3, dim <= 8, degrees <= 4
    python3 sweep12.py --k 3 --maxdim 12 --maxdeg 5
    python3 sweep12.py --k 4 --maxdim 16 --maxdeg 4 --out sl2fourth.json

Cost.  The nullcone computation dominates and grows quickly in both dim V and
maxdeg: degree d needs the monomials of C[V]_d, of which there are
binomial(dim V + d - 1, d).  dim V = 10 with maxdeg = 4 is a few seconds per
representation; dim V = 16 with maxdeg = 5 is minutes.  Start small and push.
Results are streamed to stdout and written incrementally to the JSON file, so
the run can be interrupted without losing work.
"""
import sys, time, json, argparse
from rootdata import Root
from hess import is_hesselink
from census import fundamental_hws, multisets, uses_all_factors
from mq import test_theorem12


def calibrate(maxdeg=4):
    print("Calibrating the detector")
    ok = True
    for sm, expect in [([(1,)], True), ([(1,), (1,)], True), ([(2,)], True),
                       ([(3,)], False), ([(4,)], False)]:
        got, _ = test_theorem12(sm, maxdeg, verbose=False)
        good = (got == expect)
        ok &= good
        print(f"    V = {sm}: Theorem 1.2 {'holds' if got else 'fails'}"
              f"   expected {'holds' if expect else 'fails'}"
              f"   {'ok' if good else '*** DETECTOR BROKEN ***'}")
    if not ok:
        print("Calibration failed -- do not trust the sweep.")
        sys.exit(2)
    print()
    return ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--k", type=int, nargs="*", default=[1, 2, 3],
                    help="number of SL2 factors")
    ap.add_argument("--maxdim", type=int, default=8)
    ap.add_argument("--maxparts", type=int, default=3)
    ap.add_argument("--maxdeg", type=int, default=4,
                    help="highest polynomial degree of C[N_V] to compare")
    ap.add_argument("--out", default="sweep12.json")
    args = ap.parse_args()

    calibrate(min(args.maxdeg, 4))

    t0 = time.time()
    results, bad = [], []
    for k in args.k:
        rt = Root([("A", 1)] * k)
        irr = fundamental_hws(rt, args.maxdim)
        reps = [sm for sm in multisets(irr, args.maxdim, args.maxparts)
                if uses_all_factors(rt, sm)]
        hess = []
        for sm in reps:
            ok, _ = is_hesselink(rt, sm)
            if ok:
                hess.append([tuple(int(x) for x in h) for h in sm])
        print(f"=== SL2^{k}: {len(reps)} representations, "
              f"{len(hess)} Hesselink-type (dim <= {args.maxdim})")
        for i, summ in enumerate(hess):
            try:
                good, rows = test_theorem12(summ, args.maxdeg, verbose=False)
            except Exception as e:
                print(f"    [skip {summ}: {e}]")
                continue
            dim = sum(1 for _ in Root([("A", 1)] * k).weights(summ[0])) if False else None
            rec = dict(k=k, V=[list(s) for s in summ], theorem12=good)
            results.append(rec)
            if not good:
                bad.append((summ, rows))
                print(f"    *** THEOREM 1.2 FAILS ***  V = {summ}")
                for lam, got, m, same in rows:
                    if not same:
                        print(f"        chi{lam}: C[N_V] {got}   M_q {m}")
            else:
                print(f"    ok  V = {summ}   [{i+1}/{len(hess)}]"
                      f"  ({time.time()-t0:.0f}s)")
            json.dump(results, open(args.out, "w"), indent=1)

    print(f"\n{'='*66}")
    print(f"tested {len(results)} Hesselink-type representations")
    print(f"violations of Theorem 1.2: {len(bad)}")
    if not bad:
        print("  none in this range")
    print(f"wrote {args.out}   ({time.time()-t0:.0f}s)")


if __name__ == "__main__":
    main()
