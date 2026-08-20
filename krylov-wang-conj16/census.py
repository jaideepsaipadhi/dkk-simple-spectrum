"""
census.py -- search for counterexamples to Conjecture 1.6 of Krylov-Wang
(arXiv:2608.03314).

Conjecture 1.6.  Let G be semisimple.  If V is a Hesselink-type representation
then V is cofree, and identity (2) / Theorem 1.2 hold for V.

This script enumerates representations V of a semisimple group G, decides
exactly whether V is Hesselink-type, and for every Hesselink-type V tests
identity (2)

                            dim_q C[V^{rho^vee}]^h
                     X_V = ------------------------
                                dim_q C[V]^G

as a power-series identity.  A Hesselink-type V that fails (2) refutes the
conjecture; those are printed under VIOLATIONS and written to the output JSON.

Usage
-----
    python3 census.py                      # default sweep, a few minutes
    python3 census.py --groups A1 A1^2 A2  # choose groups
    python3 census.py --maxdim 20 --maxparts 4 --order 10
    python3 census.py --groups A1^3 --maxdim 24 --out sl2cubed.json

Group names: A1, A2, ..., B2, B3, C3, D4, G2, F4, and products written with
'^' or 'x', e.g.  A1^3,  A1xA2,  A1^2xB2.

Cost.  X_V is exact and cheap.  The expensive part is the identity-(2) test,
whose cost grows with --order and with dim V; it only runs on Hesselink-type
representations, which are rare, so raising --maxdim is usually affordable.
"""
import sys, time, itertools, json, argparse
from fractions import Fraction as F
from rootdata import Root
from hess import is_hesselink, rep_weights, split, pstr
from series import test_identity2


# ------------------------------------------------------------ group parsing
def parse_group(name):
    """'A1^2xB2' -> [('A',1),('A',1),('B',2)]"""
    comps = []
    for part in name.replace("*", "x").split("x"):
        part = part.strip()
        if not part:
            continue
        if "^" in part:
            base, mult = part.split("^")
            mult = int(mult)
        else:
            base, mult = part, 1
        typ, rank = base[0].upper(), int(base[1:])
        comps += [(typ, rank)] * mult
    return comps


def fundamental_hws(rt, maxdim):
    """All nontrivial dominant weights whose irrep has dimension <= maxdim."""
    out = []
    n = rt.n
    cap = maxdim
    for coeffs in itertools.product(range(0, cap + 1), repeat=n):
        if all(c == 0 for c in coeffs):
            continue
        if sum(coeffs) > 2 * cap:
            continue
        hw = tuple(F(c) for c in coeffs)
        try:
            d = sum(rt.weights(hw).values())
        except Exception:
            continue
        if d <= maxdim:
            out.append((hw, d))
    out.sort(key=lambda z: (z[1], z[0]))
    return out


def multisets(items, maxtotal, maxparts):
    res = []
    def rec(start, cur, tot):
        if cur:
            res.append(list(cur))
        if len(cur) >= maxparts:
            return
        for i in range(start, len(items)):
            hw, d = items[i]
            if tot + d <= maxtotal:
                cur.append(hw)
                rec(i, cur, tot + d)
                cur.pop()
    rec(0, [], 0)
    return res


def uses_all_factors(rt, summands):
    """Discard V on which some simple factor acts trivially: those are really
    representations of a smaller group and are covered by that group's sweep."""
    return all(any(hw[i] != 0 for hw in summands) for i in range(rt.n))


# --------------------------------------------------------------- the sweep
def run_group(name, maxdim, maxparts, order, verbose=True):
    comps = parse_group(name)
    rt = Root(comps)
    irr = fundamental_hws(rt, maxdim)
    reps = [sm for sm in multisets(irr, maxdim, maxparts) if uses_all_factors(rt, sm)]
    if verbose:
        print(f"\n=== G = {name}   ({len(irr)} irreps of dim <= {maxdim}, "
              f"{len(reps)} representations to test)")
    hits, viol = [], []
    t0 = time.time()
    for idx, sm in enumerate(reps):
        if verbose and idx and idx % 500 == 0:
            print(f"    ... {idx}/{len(reps)}  ({time.time()-t0:.0f}s)")
        try:
            ok, X = is_hesselink(rt, sm)
        except Exception as e:
            print(f"    [skip {sm}: {e}]")
            continue
        if not ok:
            continue
        holds, lhs, rhs, xp, _ = test_identity2(rt, sm, order)
        dim = sum(sum(rt.weights(hw).values()) for hw in sm)
        rec = dict(group=name, V=[[int(x) for x in hw] for hw in sm], dim=dim,
                   XV=xp, identity2=holds)
        hits.append(rec)
        if not holds:
            viol.append(rec)
        if verbose:
            tag = "(2) ok" if holds else "*** (2) FAILS ***"
            print(f"    HESSELINK  V = {[tuple(int(x) for x in h) for h in sm]}"
                  f"  dim {dim}   X_V = {xp}   {tag}")
    if verbose:
        print(f"    -> {len(hits)} Hesselink-type, {len(viol)} violating (2)"
              f"   ({time.time()-t0:.0f}s)")
    return hits, viol


DEFAULT_GROUPS = ["A1", "A1^2", "A2", "A1^3", "B2", "G2", "A1xA2", "A3"]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--groups", nargs="*", default=DEFAULT_GROUPS)
    ap.add_argument("--maxdim", type=int, default=12,
                    help="bound on dim V (default 12)")
    ap.add_argument("--maxparts", type=int, default=3,
                    help="max number of irreducible summands (default 3)")
    ap.add_argument("--order", type=int, default=8,
                    help="truncation order in q for the identity-(2) test")
    ap.add_argument("--out", default="census.json")
    args = ap.parse_args()

    t0 = time.time()
    allhits, allviol = [], []
    for g in args.groups:
        h, v = run_group(g, args.maxdim, args.maxparts, args.order)
        allhits += h
        allviol += v

    print(f"\n{'='*70}")
    print(f"total Hesselink-type representations found: {len(allhits)}")
    print(f"VIOLATIONS of identity (2): {len(allviol)}")
    for r in allviol:
        print(f"    {r['group']}  V = {r['V']}  dim {r['dim']}  X_V = {r['XV']}")
    if not allviol:
        print("    none in this range -- Conjecture 1.6 survives here")
    json.dump(allhits, open(args.out, "w"), indent=1)
    print(f"\nwrote {args.out}   ({time.time()-t0:.0f}s)")


if __name__ == "__main__":
    main()
