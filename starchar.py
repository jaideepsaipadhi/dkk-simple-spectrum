"""
(*) WITHOUT THE CONSTRUCTION ORDER.

The conjecture was reduced (upsets.py) to

    (*)  at every colour-t socle step, X is an up-set at the chain neighbour t-1.

X at that moment is V(w') for a prefix w' of the reduced word, hence is itself
an EXTREMAL module; and the step being colour t means c_t = <alpha_t^vee, mu'>
is positive.  So (*) is implied by the order-free statement

    (*o) every extremal module M with <alpha_t^vee, lambda - dim M> >= 1
         is an up-set at t-1.

This matters: (*o) mentions no word and no interleaving, which was the sole
obstruction left in chainG.py / staircase.py.  This script tests (*o) on every
extremal module in range, and reports the up-set failures at t-1 together with
their c_t, to see whether c_t >= 1 is exactly the right side condition.
"""
import sys, time
from collections import Counter
from rootsys import RootSystem, orbit_with_words
from wordtrace import trace
from homological import setup
from gls import socle_step
from lemA import TRIV
from upsets import upset

CASES = [("D", 5, 2), ("D", 5, 3), ("D", 6, 2), ("D", 6, 3), ("D", 6, 4),
         ("D", 7, 2), ("D", 7, 3), ("D", 7, 4), ("D", 7, 5), ("D", 8, 2),
         ("D", 8, 3), ("D", 8, 4), ("D", 9, 3), ("D", 9, 4), ("D", 10, 3)]


def run(typ, n, k, cap=3000):
    rs, pr, kk = setup(typ, n, k)
    t = TRIV[typ](n)
    lam = rs.fundamental(k)
    orb = orbit_with_words(rs, lam)
    if len(orb) > cap:
        return None
    st = Counter()
    seen = set()
    for mu, word in orb.items():
        X = {}
        mucur = list(lam)
        for _, i, c, _ in trace(rs, word, lam):
            key = tuple(sorted((d, j, len(v)) for (d, j), v in X.items() if v))
            if key not in seen:
                seen.add(key)
                ct = rs.pairing(rs.simple_root(t), mucur) if False else None
                # c_t = <alpha_t^vee, mu'> with mu' = lambda - dim X (root coords)
                vv = {j: sum(len(X.get((d, j), [])) for d in range(pr.top + 1))
                      for j in rs.I}
                ct = (1 if t == kk else 0) + sum(vv[j] for j in rs.adj[t]) - 2 * vv[t]
                u = upset(pr, X, t - 1)
                st['tot'] += 1
                st['up'] += u
                st['pos_tot'] += (ct >= 1)
                st['pos_up'] += (ct >= 1 and u)
                if not u:
                    st['fail_ct=%d' % ct] += 1
            for _ in range(c):
                X = socle_step(pr, rs, X, i)
    return st


if __name__ == "__main__":
    t0 = time.time()
    S = Counter()
    for typ, n, k in CASES:
        r = run(typ, n, k)
        if r:
            S += r
            print(f"  {typ}{n} kappa={k}: up-set at t-1 {r['up']}/{r['tot']}, "
                  f"with c_t>=1 {r['pos_up']}/{r['pos_tot']}")
    print(f"\n  ALL extremal modules, up-set at t-1 : {S['up']}/{S['tot']}")
    print(f"  those with c_t >= 1  (this is (*o))  : {S['pos_up']}/{S['pos_tot']}")
    for kk, vv in sorted(S.items()):
        if kk.startswith('fail'):
            print(f"     up-set failures with {kk[5:]}: {vv}")
    ok = S['pos_up'] == S['pos_tot']
    print(f"\n({time.time()-t0:.0f}s)")
    print("(*o) HOLDS" if ok else "*** (*o) fails ***")
    sys.exit(0 if ok else 1)
