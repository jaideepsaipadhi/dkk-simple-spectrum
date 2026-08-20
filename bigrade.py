"""
THE BIGRADING IS DETERMINED BY THE PATH-LENGTH GRADING.

Pi is bigraded by the number of a-arrows and the number of a*-arrows
separately (the relations sum_a (a a* - a* a) are bihomogeneous of bidegree
(1,1)).  A torus that sees only one of the two bidegrees would in principle
induce a decomposition of V finer or coarser than the path-length grading, and
the criterion is stated in terms of path length -- so the two must be compared.

They agree, for a reason special to Dynkin diagrams: Gamma is a TREE.

  Fix an orientation of Gamma.  Since Gamma is a tree, there is a function
  f : vertices -> Z with f(head) - f(tail) = 1 for every oriented edge (build
  it by walking the tree from any root; no cycle can obstruct it).  An
  a-arrow s -> t raises f by 1 and its reverse a* lowers it by 1, so any path
  from j to i satisfies

        (#a-arrows) - (#a*-arrows)  =  f(i) - f(j),

  determined by the endpoints alone.  Together with
  (#a) + (#a*) = path length = d, this pins the bidegree:

        p = (d + f(i) - f(j))/2,      q = (d - f(i) + f(j))/2.

  So (Pi e_j)_i[d] is PURE of a single bidegree, and both p and q are affine
  functions of d once i and j are fixed.

Consequence: for any subtorus acting through the bigrading -- whether it sees
p, or q, or p+q -- the weight on (I_kappa)_i[d] is an injective function of d.
So the weight decomposition of V_i is EXACTLY the path-length decomposition,
in both directions, and the criterion "dim V_i[d] <= 1 for all i,d" is
independent of which convention is used.

This script checks that Gamma is a tree and that f exists, for every ADE
diagram in range.
"""
import sys
from rootsys import RootSystem


def potential(rs):
    """f with f(head) - f(tail) = 1 for each edge of a fixed orientation.
    Orient every edge i -> j with i < j.  Returns None if no such f exists."""
    f = {}
    root = min(rs.I)
    f[root] = 0
    front = [root]
    while front:
        nxt = []
        for x in front:
            for y in rs.adj[x]:
                val = f[x] + (1 if x < y else -1)
                if y in f:
                    if f[y] != val:
                        return None
                else:
                    f[y] = val
                    nxt.append(y)
        front = nxt
    return f if len(f) == len(rs.I) else None


def is_tree(rs):
    edges = sum(len(rs.adj[i]) for i in rs.I) // 2
    return edges == len(rs.I) - 1


if __name__ == "__main__":
    ok = True
    print("Gamma is a tree, so the bidegree is determined by (vertex, degree)\n")
    for typ, ns in [("A", range(2, 9)), ("D", range(4, 10)), ("E", (6, 7, 8))]:
        for n in ns:
            rs = RootSystem(typ, n)
            tree = is_tree(rs)
            f = potential(rs)
            good = tree and f is not None
            ok &= good
            print(f"  {typ}{n}: tree={tree}  potential f={None if f is None else [f[i] for i in rs.I]}"
                  f"   {'OK' if good else 'FAIL'}")
    print("\n" + ("BIDEGREE IS A FUNCTION OF (VERTEX, PATH-LENGTH DEGREE)"
                  if ok else "*** failed ***"))
    sys.exit(0 if ok else 1)
