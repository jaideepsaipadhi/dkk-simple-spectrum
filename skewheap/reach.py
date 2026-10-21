"""Can Proctor's hook recursion (Def 2.6) reach the TRUE hook multiset at all?

For a hook length poset P the skew formula at F = empty reduces to
    series(P) = prod_u 1/(1 - q^{h(u)}),
so any hook function the formula could use must have the true multiset.  Outside the
d-complete class an element may top several d_k-intervals, so Def 2.6 is a FAMILY of
recursions, one per choice.  This script asks, for every hook-length non-d-complete heap,
whether ANY choice in that family produces the true multiset.

If none ever does, the Naruse-Okada formula fails on these posets at F = empty already --
before excited diagrams enter at all.
usage: python3 reach.py [MAXLEN=10] [PMAX=11] [TYPES=D4,D5]
"""
import sys, itertools, collections
sys.path.insert(0, '.')
from nsl import cartan, heap_down, is_reduced
from dcomp2 import dcomplete
from novel import dual, relabel_natural
from no26 import leq_matrix, dk_intervals
from filterhook import hooks_of
from skew import to_proctor
from choice import hooks_from_choice

def main():
    ML   = int(sys.argv[1]) if len(sys.argv) > 1 else 10
    PMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 11
    TY   = (sys.argv[3] if len(sys.argv) > 3 else "D4,D5").split(",")
    stat = collections.Counter(); ex = []
    for t in TY:
        typ, rk = t[0], int(t[1:]); C = cartan(typ, rk)
        out = []
        def rec(x):
            if x and is_reduced(C, x): out.append(tuple(x))
            if len(x) == ML: return
            for i in range(rk): rec(x + [i])
        rec([])
        for w in out:
            down = heap_down(C, list(w)); n = len(down)
            if n < 4 or n > PMAX: continue
            dn = relabel_natural(down); th = hooks_of(dn)
            if th is None: continue                       # hook length posets only
            dc = (len(dcomplete(relabel_natural(dual(down)))) == 0)
            dP, col, order = to_proctor(down, list(w)); leq = leq_matrix(dP)
            tops = collections.defaultdict(list)
            for (v, u, x, y) in dk_intervals(dP, leq): tops[u].append((v, x, y))
            keys = sorted(tops)
            nc = 1
            for u in keys: nc *= (len(tops[u]) + 1)
            if nc > 20000: stat[(dc, 'too many choices')] += 1; continue
            reach = False
            for pick in itertools.product(*[[None]+list(range(len(tops[u]))) for u in keys]):
                ch = {u: tops[u][p] for u, p in zip(keys, pick) if p is not None}
                if sorted(hooks_from_choice(dP, leq, ch)) == th: reach = True; break
            stat[('d-complete' if dc else 'not d-complete',
                  f'true multiset reachable: {reach}')] += 1
            if dc and not reach and len(ex) < 3: ex.append(('D-COMPLETE UNREACHABLE', t, w, th))
            if (not dc) and reach and len(ex) < 3: ex.append(('NON-DC REACHABLE', t, w, th))
    print(f"types {TY}, words to length {ML}, |P| <= {PMAX}; hook length posets only\n")
    for k in sorted((k for k in stat if isinstance(k, tuple)), key=str):
        print(f"   {k[0]:16} {k[1]:30} {stat[k]}")
    for a in ex: print(f"\n   {a}")

if __name__ == "__main__":
    main()
