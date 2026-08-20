"""
Candidate level-assignment rule (the missing "generalized heap").

Process the reduced word right-to-left.  At the step with letter i and
multiplicity c = <alpha_i^vee, mu>:

    candidates = { l + 1 : l a level already occupied by a NEIGHBOUR colour }
    place the c lowest candidates not already occupied by colour i.

(The first step, which has no placed neighbours, goes to level 0.)

Tested against the three posets printed in DKK Sec 5.4.2.
"""

from rootsys import RootSystem
from wordtrace import trace


def build_levels(rs, word, lam, defer=True):
    """Right-to-left sweep.

    At a step with letter i and multiplicity c: place ONE element at the
    lowest free candidate (neighbour level + 1).  The remaining c-1 elements
    are DEFERRED to the end of the sweep, then placed above the final
    neighbour chain.  Deferring matters: at the moment of a c=2 step the
    neighbour chain has not finished growing, and the extra element sits on
    top of its final height, not its height at that instant.
    """
    L = {i: [] for i in rs.I}
    pending = []

    def candidates(i):
        nb = sorted({l for j in rs.adj[i] for l in L[j]})
        return [l + 1 for l in nb] if nb else [0]

    def place_one(i):
        for x in candidates(i):
            if x not in L[i]:
                L[i].append(x)
                L[i].sort()
                return
        L[i].append((max(L[i]) + 2) if L[i] else 0)
        L[i].sort()

    for _, i, c, _ in trace(rs, word, lam):
        if c >= 1:
            place_one(i)
        if defer:
            pending += [i] * (c - 1)
        else:
            for _ in range(c - 1):
                place_one(i)

    for i in pending:                       # now the chains are final
        nb = [l for j in rs.adj[i] for l in L[j]]
        x = (max(nb) + 1) if nb else 0
        while x in L[i]:
            x += 1
        L[i].append(x)
        L[i].sort()

    return {i: L[i] for i in rs.I if L[i]}


def levels_to_rows(L):
    top = max(max(v) for v in L.values())
    return [sorted(i for i, ls in L.items() if d in ls) for d in range(top + 1)]


def check(label, typ, n, k, word, paper_rows):
    rs = RootSystem(typ, n)
    L = build_levels(rs, word, rs.fundamental(k))
    rows = levels_to_rows(L)
    ok = rows == paper_rows
    print(f"{label}  {rs.name()}  w={''.join('s'+str(i) for i in word)}")
    print(f"    computed levels : {dict(sorted(L.items()))}")
    print(f"    computed rows   : {rows}")
    print(f"    paper rows      : {paper_rows}")
    print(f"    {'MATCH' if ok else 'MISMATCH'}\n")
    return ok


if __name__ == "__main__":
    print("Testing the candidate level rule against DKK's three posets\n")
    r = []
    r.append(check("Ex 5.10", "D", 4, 2, (1, 2, 3, 4, 2),
                   [[2], [1, 3, 4], [2], [1]]))
    r.append(check("Ex 5.11", "D", 5, 2, (2, 1, 3, 4, 5, 3, 2),
                   [[2], [1, 3], [2, 4, 5], [3], [2]]))
    r.append(check("Ex 5.12", "D", 5, 2, (2, 1, 2, 3, 4, 5, 3, 2),
                   [[2], [1, 3], [2, 4, 5], [3], [2], [1]]))
    print(f"{sum(r)}/3 reproduced")
