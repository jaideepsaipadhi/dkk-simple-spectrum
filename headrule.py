"""
LEMMA H' (the head of an iterated socle is the set of "heap-maximal" colours).

For X = soc_sigma(I) built from a colour sequence sigma:

        Hom(X, S_l) != 0   <=>   l occurs in sigma and NO neighbour of l
                                 occurs in sigma after l's last occurrence.

Only the "<=" direction (contrapositive) is needed downstream:

  (*)  if some neighbour of i occurs after i's last occurrence in sigma,
       then top_i(X) = 0.

Together with the elementary c-trace fact that between two consecutive
colour-i steps some neighbouring colour must intervene (only reflections at
i and its neighbours change <alpha_i^vee, mu>, and the value flips from c to
-c at a colour-i step, so it cannot return to >= 1 unaided), (*) gives
top_i(X) = 0 at every step, hence -- with the Crawley-Boevey identity --
dim soc_i(I/X) = c.  That is Lemma S.

This script tests H' in both directions, and also the multiplicity.
"""
import sys, time
from collections import Counter
from rootsys import RootSystem, orbit_with_words
from wordtrace import trace
from proj2 import Proj2
from gls import socle_step
from homological import setup, top_i


def predicted_head(rs, seq):
    """{l : l in seq and no neighbour of l after l's last occurrence}"""
    lastpos = {}
    for p, l in enumerate(seq):
        lastpos[l] = p
    out = set()
    for l, p in lastpos.items():
        if not any(x in rs.adj[l] for x in seq[p + 1:]):
            out.add(l)
    return out


def run(typ, n, k, cap=300):
    rs, pr, kk = setup(typ, n, k)
    lam = rs.fundamental(k)
    orb = orbit_with_words(rs, lam)
    if len(orb) > cap:
        return None
    tot = ok = 0
    star_tot = star_ok = 0
    bad = Counter()
    for mu, word in orb.items():
        if not word:
            continue
        X, seq = {}, []
        for _, i, c, _ in trace(rs, word, lam):
            if seq:
                heads = {l for l in rs.I if top_i(pr, rs, X, l)}
                pred = predicted_head(rs, seq)
                tot += 1
                ok += (heads == pred)
                if heads != pred:
                    bad[(typ, n, k, tuple(sorted(heads)), tuple(sorted(pred)))] += 1
                # the direction actually used downstream
                star_tot += 1
                star_ok += (i not in heads)
            for _ in range(c):
                X = socle_step(pr, rs, X, i)
                seq.append(i)
    return tot, ok, star_tot, star_ok, bad


CASES = [("A", 4, 2), ("A", 5, 3), ("D", 4, 2), ("D", 5, 2), ("D", 5, 3),
         ("E", 6, 2), ("E", 6, 3), ("E", 7, 1)]

if __name__ == "__main__":
    t0 = time.time()
    T = O = ST = SO = 0
    B = Counter()
    for typ, n, k in CASES:
        r = run(typ, n, k)
        if r is None:
            continue
        t, o, st, so, b = r
        T += t; O += o; ST += st; SO += so; B += b
        print(f"  {typ}{n} om{k}: H' exact {o}/{t};  (*) top_i(X)=0 {so}/{st}")
    print(f"\ntotals: H' {O}/{T},  (*) {SO}/{ST}   ({time.time()-t0:.0f}s)")
    if B:
        print("H' deviations (got head, predicted head):")
        for kk_, vv in B.most_common(8):
            print("   ", kk_, vv)
    sys.exit(0 if SO == ST else 1)
