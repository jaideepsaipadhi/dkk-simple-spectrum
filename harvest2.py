"""
Harvest v2 -- corrected indexing set.

mu = w.lambda depends only on the coset w W_lambda, so the right index set is
W^lambda, the minimal-length coset representatives, equivalently the orbit
W.lambda itself. Enumerating all of W (harvest.py) produced spurious c=0
steps from non-minimal representatives.

For each mu in the orbit we take a shortest word w with w.lambda = mu and
record the c-trace.
"""

from collections import Counter
from rootsys import RootSystem, orbit_with_words
from wordtrace import trace, v_from_trace
from weylprops import is_fully_commutative


def study(typ, n, k, do_fc=True, show=6):
    rs = RootSystem(typ, n)
    lam = rs.fundamental(k)
    orb = orbit_with_words(rs, lam)
    print(f"=== {rs.name()}   lambda = omega_{k}   |orbit| = {len(orb)}")

    maxc = Counter()
    zero = 0
    defects = []            # elements with some c >= 2
    for mu, word in sorted(orb.items(), key=lambda t: len(t[1])):
        if not word:
            continue
        tr = trace(rs, word, lam)
        cs = [c for _, _, c, _ in tr]
        if min(cs) == 0:
            zero += 1
        m = max(cs)
        maxc[m] += 1
        if m >= 2:
            v = v_from_trace(rs, tr)
            defects.append((word, cs, v, sum(1 for c in cs if c >= 2)))

    print(f"  steps with c=0 (should be 0 now): {zero}")
    print(f"  c_max distribution: "
          + ", ".join(f"{m}:{c}" for m, c in sorted(maxc.items())))
    if defects:
        ndef = Counter(d[3] for d in defects)
        print(f"  number of defect steps (c>=2) per element: "
              + ", ".join(f"{k2} defect(s):{v2}" for k2, v2 in sorted(ndef.items())))
        print(f"  sample defective elements:")
        for word, cs, v, nd in defects[:show]:
            s = ''.join('s' + str(i) for i in word)
            fc = is_fully_commutative(rs, word) if do_fc and len(word) <= 12 else "-"
            print(f"    w={s:<28} c={cs}  v={v}  FC={fc}")
    print()
    return maxc


if __name__ == "__main__":
    print("Harvest v2: minimal coset representatives (orbit of lambda)\n")
    tot = Counter()
    for typ, n, k in [("D", 4, 2), ("D", 5, 2), ("D", 6, 2), ("D", 6, 3),
                      ("A", 4, 2), ("D", 4, 1), ("E", 6, 2), ("E", 6, 1)]:
        try:
            tot += study(typ, n, k, do_fc=(n <= 6))
        except Exception as e:
            print(f"  ({typ}{n}, om{k}) failed: {e}\n")
    print("Aggregate c_max distribution over all runs:",
          dict(sorted(tot.items())))
    print("\nc_max never exceeds 2?", max(tot) <= 2)
