"""
Type E: the graded multiplicity profile of Pi e_k, for every fundamental
weight of E6, E7, E8.  This is a FINITE table -- unlike type D there is no
rank to induct on -- so establishing the type-E analogues of Lemmas P and A
is a finite verification, not an induction.

Reports, for each (type, k):
  * the profile of the trivalent node t
  * the max multiplicity anywhere
  * WINDOW CONTAINMENT (Lemma A): is  {d : mult_i(d) >= 2}  contained in
    {d : mult_t(d) >= 2}  for every vertex i?
"""
from rootsys import RootSystem
from proj2 import Proj2

TRIV = 4   # trivalent node in the E-labelling of rootsys.py


def profile(pr, rs, i):
    return {d: pr.dim[(d, i)] for (d, j) in pr.dim if j == i}


def report(typ, n):
    rs = RootSystem(typ, n)
    print(f"=== {rs.name()}  (trivalent node t = {TRIV}) ===")
    allok = True
    for k in rs.I:
        pr = Proj2(rs, k, dmax=40)
        pt = profile(pr, rs, TRIV)
        big_t = {d for d, m in pt.items() if m >= 2}
        mx = max(pr.dim.values())
        bad = []
        for i in rs.I:
            bi = {d for d, m in profile(pr, rs, i).items() if m >= 2}
            if not bi <= big_t:
                bad.append((i, sorted(bi - big_t)))
        allok &= not bad
        seq = [pt[d] for d in sorted(pt)]
        degs = sorted(pt)
        print(f"  k={k}: (Pi e_k)_t = {seq}  at degrees {degs}"
              f"   maxmult={mx}   Lemma A: {'OK' if not bad else bad}")
    return allok


if __name__ == "__main__":
    ok = True
    for typ, n in [("E", 6), ("E", 7), ("E", 8)]:
        ok &= report(typ, n)
        print()
    print("WINDOW CONTAINMENT HOLDS IN ALL OF TYPE E" if ok
          else "*** containment fails somewhere ***")
