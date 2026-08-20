"""
The generalized heap: levels AND covering relations.

Sweep the reduced word right-to-left.  At the step with letter i and
multiplicity c = <alpha_i^vee, mu>:

  * place ONE element of colour i at the lowest free candidate level, where
    candidates = { l+1 : l a level occupied by an ALREADY PLACED neighbour }
    (level 0 if no neighbour is placed yet);
  * DEFER the remaining c-1 elements to the end of the sweep, then place them
    one level above the FINAL neighbour chain.

EDGE RULE: a new element at level x covers exactly those ALREADY EXISTING
elements of neighbouring colours at level x-1.  Elements created later never
become covered by elements created earlier -- this is what deletes the edges
that a naive "adjacent colour, consecutive level" reading would insert, and
it is what makes the moment-map relations consistent.

Verified against the Hasse diagrams printed on p.82-83 of arXiv:2608.16746.
"""

from rootsys import RootSystem
from wordtrace import trace


class Heap:
    def __init__(self, rs):
        self.rs = rs
        self.elts = []            # list of (colour, level), index = identity
        self.cov = []             # cov[a] = set of b with b <| a  (a covers b)

    def levels_of(self, colour):
        return [l for (c, l) in self.elts if c == colour]

    def add(self, colour, level):
        idx = len(self.elts)
        below = {n for n, (c, l) in enumerate(self.elts)
                 if l == level - 1 and c in self.rs.adj[colour]}
        self.elts.append((colour, level))
        self.cov.append(below)
        return idx

    def rows(self):
        top = max(l for _, l in self.elts)
        return [sorted(c for c, l in self.elts if l == d) for d in range(top + 1)]

    def edges(self):
        return sorted((b, a) for a, bs in enumerate(self.cov) for b in bs)

    def below_closure(self, a):
        seen, stack = set(), [a]
        while stack:
            x = stack.pop()
            for b in self.cov[x]:
                if b not in seen:
                    seen.add(b); stack.append(b)
        return seen


def build_heap(rs, word, lam):
    H = Heap(rs)

    def candidates(colour):
        nb = sorted({l for (c, l) in H.elts if c in rs.adj[colour]})
        return [l + 1 for l in nb] if nb else [0]

    def place(colour):
        used = H.levels_of(colour)
        for x in candidates(colour):
            if x not in used:
                return H.add(colour, x)
        return H.add(colour, (max(used) + 2) if used else 0)

    pending = []
    for _, i, c, _ in trace(rs, word, lam):
        if c >= 1:
            place(i)
        pending += [i] * (c - 1)

    for i in pending:
        nb = [l for (c, l) in H.elts if c in rs.adj[i]]
        x = (max(nb) + 1) if nb else 0
        used = H.levels_of(i)
        while x in used:
            x += 1
        H.add(i, x)
    return H


# --------------------------------------------------- moment-map consistency

def relation_defects(H):
    """Same-coloured pairs x < z with lev(z)=lev(x)+2 and exactly ONE
    intermediate: the preprojective relation then forces a zero map, so the
    configuration is NOT realizable."""
    bad = []
    n = len(H.elts)
    for a in range(n):
        ca, la = H.elts[a]
        for b in range(n):
            cb, lb = H.elts[b]
            if cb != ca or lb != la + 2:
                continue
            mids = [m for m in H.cov[b] if a in H.cov[m]]
            if len(mids) == 1:
                bad.append((H.elts[a], H.elts[b], len(mids)))
    return bad


def simple_spectrum(H):
    """Distinct (colour, level) pairs -- i.e. each V_i has multiplicity-free
    torus spectrum."""
    seen = set()
    for e in H.elts:
        if e in seen:
            return False
        seen.add(e)
    return True


# ------------------------------------------------------------------- tests

PAPER = {
    "Ex 5.10": (("D", 4, 2), (1, 2, 3, 4, 2),
                [[2], [1, 3, 4], [2], [1]],
                # edges as (lower colour@level -> upper colour@level)
                {((2, 0), (3, 1)), ((2, 0), (1, 1)), ((2, 0), (4, 1)),
                 ((3, 1), (2, 2)), ((4, 1), (2, 2)), ((2, 2), (1, 3))}),
    "Ex 5.11": (("D", 5, 2), (2, 1, 3, 4, 5, 3, 2),
                [[2], [1, 3], [2, 4, 5], [3], [2]],
                {((2, 0), (3, 1)), ((2, 0), (1, 1)),
                 ((3, 1), (5, 2)), ((3, 1), (4, 2)), ((3, 1), (2, 2)),
                 ((1, 1), (2, 2)),
                 ((5, 2), (3, 3)), ((4, 2), (3, 3)), ((3, 3), (2, 4))}),
    "Ex 5.12": (("D", 5, 2), (2, 1, 2, 3, 4, 5, 3, 2),
                [[2], [1, 3], [2, 4, 5], [3], [2], [1]],
                {((2, 0), (3, 1)), ((2, 0), (1, 1)),
                 ((3, 1), (5, 2)), ((3, 1), (4, 2)), ((3, 1), (2, 2)),
                 ((1, 1), (2, 2)),
                 ((5, 2), (3, 3)), ((4, 2), (3, 3)),
                 ((3, 3), (2, 4)), ((2, 4), (1, 5))}),
}


if __name__ == "__main__":
    print("Generalized heap vs the Hasse diagrams on p.82-83\n")
    allok = True
    for label, ((typ, n, k), word, rows, edges) in PAPER.items():
        rs = RootSystem(typ, n)
        H = build_heap(rs, word, rs.fundamental(k))
        got_rows = H.rows()
        got_edges = {(H.elts[b], H.elts[a]) for (b, a) in H.edges()}
        okr = got_rows == rows
        oke = got_edges == edges
        allok &= okr and oke
        print(f"{label}  {rs.name()}")
        print(f"   rows  {'OK ' if okr else 'BAD'} {got_rows}")
        print(f"   edges {'OK ' if oke else 'BAD'} ({len(got_edges)} edges)")
        if not oke:
            print(f"      missing : {sorted(edges - got_edges)}")
            print(f"      spurious: {sorted(got_edges - edges)}")
        print(f"   moment-map defects: {relation_defects(H)}")
        print(f"   simple spectrum   : {simple_spectrum(H)}")
        print()
    print("ALL THREE HEAPS REPRODUCED (rows + edges)" if allok else "MISMATCH")
