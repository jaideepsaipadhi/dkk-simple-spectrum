"""
Corrected base rule, checked against ground truth.

STANDARD RULE (provably right when every c = 1): sweeping the word
right-to-left, an element's level is
      1 + max( levels of all already-placed elements whose colour does NOT
               commute with it -- i.e. equal or adjacent colour )
or 0 if none.  Note SAME colour counts: that is what my old "lowest free
candidate" rule was missing, and it is why every failure was low by exactly 2.
"""
from rootsys import RootSystem, orbit_with_words
from wordtrace import trace
from trueheap import true_heap, is_lambda_min


def build_levels_std(rs, word, lam):
    L = {i: [] for i in rs.I}
    placed = []                              # (colour, level)
    pending = []
    def lev_for(i):
        rel = [l for (c, l) in placed if c == i or c in rs.adj[i]]
        return (max(rel) + 1) if rel else 0
    for _, i, c, _ in trace(rs, word, lam):
        x = lev_for(i); L[i].append(x); placed.append((i, x))
        pending += [i] * (c - 1)
    for i in pending:
        x = lev_for(i)
        while x in L[i]:
            x += 1
        L[i].append(x); placed.append((i, x))
    return {i: sorted(v) for i, v in L.items() if v}


ok = tot = 0
bad = []
for typ, n, k in [("A",4,2),("A",5,3),("D",4,1),("D",5,1),("D",6,1),
                  ("E",6,1),("D",4,2),("D",5,2),("D",6,2),("E",6,2),
                  ("D",7,3),("E",7,1)]:
    rs = RootSystem(typ, n); lam = rs.fundamental(k)
    o = t = 0
    for mu, word in orbit_with_words(rs, lam).items():
        if not word or not is_lambda_min(rs, word, lam):
            continue
        t += 1
        if build_levels_std(rs, word, lam) == true_heap(rs, word):
            o += 1
        elif len(bad) < 3:
            bad.append((rs.name(), word))
    ok += o; tot += t
    print(f"{rs.name()} om{k}: {o}/{t}")
print(f"\nTOTAL lambda-minuscule: {ok}/{tot}", "ALL CORRECT" if ok == tot else f"  failures: {bad}")

# now the three published (non-lambda-minuscule) examples
print("\nSame rule on DKK's three examples:")
for label, (typ, n, k), word, truth in [
    ("Ex 5.10", ("D",4,2), (1,2,3,4,2), {1:[1,3],2:[0,2],3:[1],4:[1]}),
    ("Ex 5.11", ("D",5,2), (2,1,3,4,5,3,2), {1:[1],2:[0,2,4],3:[1,3],4:[2],5:[2]}),
    ("Ex 5.12", ("D",5,2), (2,1,2,3,4,5,3,2), {1:[1,5],2:[0,2,4],3:[1,3],4:[2],5:[2]}),
]:
    rs = RootSystem(typ, n)
    got = build_levels_std(rs, word, rs.fundamental(k))
    print(f"  {label}: {'OK' if got==truth else 'WRONG'}   got {dict(sorted(got.items()))}")
