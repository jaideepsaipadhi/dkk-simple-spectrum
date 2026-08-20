"""
Generalized heap, v2 -- edges determined by MOMENT-MAP CONSISTENCY.

Levels: as before (deferred-placement sweep), verified 3/3.

Edges: start from ALL candidate covers (neighbouring colour, consecutive
level).  A same-coloured pair x < z with lev(z)=lev(x)+2 and exactly ONE
intermediate y forces a_{xy} a_{yz} = 0, which is impossible.  Repair by
deleting the LOWER edge x <| y.  Iterate to a fixed point.
"""
from rootsys import RootSystem
from levels import build_levels
from heap import PAPER


def build_edges(rs, L):
    elts = sorted([(i, l) for i, ls in L.items() for l in ls],
                  key=lambda e: (e[1], e[0]))
    E = {(a, b) for a in elts for b in elts
         if b[1] == a[1] + 1 and b[0] in rs.adj[a[0]]}
    while True:
        forced = None
        for x in elts:
            for z in elts:
                if z[0] != x[0] or z[1] != x[1] + 2:
                    continue
                mids = [y for y in elts
                        if (x, y) in E and (y, z) in E]
                if len(mids) == 1:
                    forced = (x, mids[0])
                    break
            if forced:
                break
        if not forced:
            return elts, E
        E.discard(forced)


if __name__ == "__main__":
    print("Edges from moment-map consistency\n")
    ok = True
    for label, ((typ, n, k), word, rows, paper_edges) in PAPER.items():
        rs = RootSystem(typ, n)
        L = build_levels(rs, word, rs.fundamental(k))
        elts, E = build_edges(rs, L)
        match = E == paper_edges
        ok &= match
        print(f"{label}  {rs.name()}   edges {'OK' if match else 'BAD'} "
              f"({len(E)} vs paper {len(paper_edges)})")
        if not match:
            print(f"   missing : {sorted(paper_edges - E)}")
            print(f"   spurious: {sorted(E - paper_edges)}")
    print("\nALL THREE HASSE DIAGRAMS REPRODUCED" if ok else "\nMISMATCH")
