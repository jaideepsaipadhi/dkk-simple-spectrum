"""
GROUND TRUTH for lambda-minuscule w: the honest heap H(w) of Definition 5.2.

Elements = letters of a reduced word; a < b iff a comes LATER in the word and
the colours do not commute (adjacent or equal).  Level = rank function.
For lambda-minuscule w this is exactly the poset DKK use, so it is an
independent check on my constructor -- no fitting, no guessing.

Also: the CORRECT reading of Dranowski et al. Prop 2.7 / 3.8.  The relation
concerns CONSECUTIVE elements of a colour fibre (which need NOT be two levels
apart), and between them lie exactly two neighbours z1, z2:
   pi(z1) != pi(z2)  -> diamond, signs cancel  (-x + x = 0)
   pi(z1) == pi(z2)  -> single runner, composite vanishes automatically
My earlier "distance-2 with >=2 intermediates" test was wrong on both counts.
"""
from rootsys import RootSystem, orbit_with_words
from wordtrace import trace
from levels import build_levels


def true_heap(rs, word):
    """word is left-to-right; later position = lower in the heap."""
    n = len(word)
    idx = list(range(n))                      # 0..n-1, position order
    def noncommute(a, b):
        return word[a] == word[b] or word[b] in rs.adj[word[a]]
    # a < b  iff  a later in word (bigger index) and non-commuting
    less = {a: set() for a in idx}
    for a in idx:
        for b in idx:
            if a > b and noncommute(a, b):
                less[b].add(a)                # a below b
    # rank by longest chain below
    lev = {}
    for b in sorted(idx, reverse=True):       # from lowest (largest index) up
        below = [lev[a] for a in less[b] if a in lev]
        lev[b] = (max(below) + 1) if below else 0
    L = {}
    for p in idx:
        L.setdefault(word[p], []).append(lev[p])
    return {i: sorted(v) for i, v in L.items()}


def is_lambda_min(rs, word, lam):
    return all(c == 1 for _, _, c, _ in trace(rs, word, lam))


print("Constructor vs TRUE heap on lambda-minuscule elements")
print("(ground truth known independently -- no fitting)\n")
grand_ok = grand_tot = 0
for typ, n, k in [("A",4,2),("A",5,3),("D",4,1),("D",5,1),("D",6,1),
                  ("E",6,1),("D",4,2),("D",5,2),("D",6,2),("E",6,2)]:
    rs = RootSystem(typ, n); lam = rs.fundamental(k)
    orb = orbit_with_words(rs, lam)
    ok = tot = 0
    firstbad = None
    for mu, word in sorted(orb.items(), key=lambda t: len(t[1])):
        if not word or not is_lambda_min(rs, word, lam):
            continue
        tot += 1
        mine = build_levels(rs, word, lam)
        truth = true_heap(rs, word)
        if mine == truth:
            ok += 1
        elif firstbad is None:
            firstbad = (word, mine, truth)
    grand_ok += ok; grand_tot += tot
    tag = "ALL OK" if ok == tot else f"{tot-ok} WRONG"
    print(f"{rs.name()} om{k}: {ok}/{tot} lambda-minuscule elements match   {tag}")
    if firstbad:
        w, m, t = firstbad
        print(f"     w={''.join('s'+str(i) for i in w)}")
        print(f"       mine  = {dict(sorted(m.items()))}")
        print(f"       truth = {dict(sorted(t.items()))}")
print(f"\nTOTAL: {grand_ok}/{grand_tot}")
