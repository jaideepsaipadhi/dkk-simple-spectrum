"""
Independent validity check: is the constructed colored poset actually
realizable as a graded quiver representation satisfying the moment map?

Prop 5.4 model: V_i has basis {e_x : pi(x)=i}; the quiver maps send
e_x -> a_{xy} e_y whenever y covers x (all a_{xy} nonzero); both Y and Z
raise the level by 1, so the preprojective relation at vertex i is a map
V_i[d] -> V_i[d+2] and reads, for each x of colour i and each z of colour i
with lev(z) = lev(x)+2:

    sum over intermediate y (pi(y)=j ~ i, x < y < z)  eps(i,j) a_{xy} a_{yz} = 0

KEY COMBINATORIAL CONSEQUENCE: if such a pair (x,z) has exactly ONE
intermediate y, the equation forces a_{xy} a_{yz} = 0, contradicting
nonvanishing.  So

    every same-coloured pair at level distance 2 must have >= 2 intermediates

is a necessary condition for realizability -- checkable with no linear algebra.
"""

from rootsys import RootSystem, orbit_with_words
from levels import build_levels


def poset_elements(L):
    """[(colour, level)] for the constructed poset."""
    return [(i, l) for i, ls in L.items() for l in ls]


def intermediates(rs, L, i, la, lb):
    """Elements of colour adjacent to i sitting at level la+1, between."""
    return [(j, la + 1) for j in rs.adj[i] if (la + 1) in L.get(j, [])]


def bad_pairs(rs, L):
    """Same-coloured pairs at level distance 2 with fewer than 2 intermediates."""
    out = []
    for i, ls in L.items():
        for a in ls:
            if a + 2 in ls:
                mids = intermediates(rs, L, i, a, a + 2)
                if len(mids) < 2:
                    out.append(((i, a), (i, a + 2), len(mids)))
    return out


def check(label, typ, n, k, word):
    rs = RootSystem(typ, n)
    L = build_levels(rs, word, rs.fundamental(k))
    bp = bad_pairs(rs, L)
    print(f"{label}  {rs.name()}  levels={dict(sorted(L.items()))}")
    if bp:
        for (x, z, m) in bp:
            print(f"    UNREALIZABLE: colour {x[0]} at levels {x[1]},{z[1]} "
                  f"has {m} intermediate(s) -- relation forces a zero map")
    else:
        print("    all same-colour distance-2 pairs have >=2 intermediates")
    print()
    return not bp


if __name__ == "__main__":
    print("Realizability probe on DKK's three examples\n")
    check("Ex 5.10", "D", 4, 2, (1, 2, 3, 4, 2))
    check("Ex 5.11", "D", 5, 2, (2, 1, 3, 4, 5, 3, 2))
    check("Ex 5.12", "D", 5, 2, (2, 1, 2, 3, 4, 5, 3, 2))

    print("Sweep over orbits: how many elements fail the necessary condition\n")
    for typ, n, k in [("D", 4, 2), ("D", 5, 2), ("D", 6, 2), ("D", 6, 3),
                      ("E", 6, 2)]:
        rs = RootSystem(typ, n)
        lam = rs.fundamental(k)
        orb = orbit_with_words(rs, lam)
        nbad = 0
        tot = 0
        for mu, word in orb.items():
            if not word:
                continue
            tot += 1
            if bad_pairs(rs, build_levels(rs, word, lam)):
                nbad += 1
        print(f"  {rs.name()} om{k}: {nbad}/{tot} constructed posets fail "
              f"the >=2-intermediates test")
