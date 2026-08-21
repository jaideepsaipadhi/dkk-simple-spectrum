"""
KRYLOV'S FOLLOW-UP QUESTION.

When w has the simple spectrum property, is the reverse-plane-partition
generating function of the resulting poset P(w) the hook-product

              1
    ------------------------------- ,      H(mu) = { <alpha^vee, mu> :
    prod_{h in H(mu)} ( 1 - z^h )                    alpha a NEGATIVE root,
                                                     <alpha^vee, mu> > 0 } ?

(a Peterson-Proctor style formula; cf. eq. (6) of Krylov's other paper).

RPP generating function.  An order-preserving map f : P -> Z_{>=0} is the same
as a decreasing multichain of order ideals I_1 >= I_2 >= ... , via
f = sum_k 1_{I_k}; so

    RPP_P(z) = sum over multichains of nonempty order ideals of z^{sum |I_j|},

which we evaluate exactly as a truncated power series by dynamic programming
over the ideal lattice.  This is convention-free: RPP_P = RPP_{P^op}.

The poset P(w) is the generalized heap of heap.py, whose covering relations
were checked against the Hasse diagrams printed on pp. 82-83 of DKK.
"""
import sys, time
from itertools import combinations
from rootsys import RootSystem, orbit_with_words
from heap import build_heap, simple_spectrum

N = 14                                    # series truncation


def mul(a, b):
    c = [0] * N
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                if y and i + j < N:
                    c[i + j] += x * y
    return c


def add(a, b):
    return [x + y for x, y in zip(a, b)]


def roots(rs):
    """All roots, in omega coordinates, via the Weyl orbit of the simples."""
    seen, frontier = set(), []
    for i in rs.I:
        a = rs.simple_root(i)
        seen.add(a)
        frontier.append(a)
    while frontier:
        new = []
        for a in frontier:
            for i in rs.I:
                b = rs.act(i, a)
                if b not in seen:
                    seen.add(b)
                    new.append(b)
        frontier = new
    return seen


def hooks(rs, mu):
    """H(mu) = multiset of <alpha^vee, mu> over NEGATIVE roots with value > 0."""
    out = []
    for a in roots(rs):
        v = rs.root_coords(a)
        if all(x <= 0 for x in v):                    # negative root
            p = rs.pairing(a, mu)
            if p > 0:
                out.append(int(p))
    return sorted(out)


def hook_product(H):
    """prod 1/(1 - z^h) as a truncated series."""
    f = [0] * N
    f[0] = 1
    for h in H:
        g = [0] * N
        for d in range(0, N, h):
            g[d] = 1
        f = mul(f, g)
    return f


def ideals(cov_down):
    """All order ideals (downward-closed sets) as frozensets of indices."""
    n = len(cov_down)
    res = set()
    for r in range(n + 1):
        for S in combinations(range(n), r):
            S = frozenset(S)
            if all(cov_down[x] <= S for x in S):
                res.add(S)
    return sorted(res, key=len)


def geom(m):
    """1/(1 - z^m) as a truncated series."""
    g = [0] * N
    for d in range(0, N, m):
        g[d] = 1
    return g


def rpp_series(cov_down):
    """RPP generating function, truncated to degree < N.

    An order-preserving f : P -> Z_{>=0} is a decreasing MULTIchain of nonempty
    order ideals I_1 >= I_2 >= ... (repeats allowed), f = sum_k 1_{I_k}.  With
    g[I] = the sum over multichains whose first term is I,

        g[I] = z^{|I|} ( 1 + g[I] + sum_{J < I} g[J] ),

    i.e.  g[I] = z^{|I|}/(1 - z^{|I|}) * ( 1 + sum_{J strictly inside I} g[J] ).
    """
    ids = [I for I in ideals(cov_down) if I]
    g = {}
    tot = [0] * N
    tot[0] = 1
    for I in ids:                                  # increasing size
        s = [0] * N
        s[0] = 1
        for J in ids:
            if len(J) < len(I) and J <= I:
                s = add(s, g[J])
        m = len(I)
        zI = [0] * N
        if m < N:
            zI[m] = 1
        g[I] = mul(mul(zI, geom(m)), s)
        tot = add(tot, g[I])
    return tot


def down_closure(H):
    """cov[a] is the set of elements a covers; return full down-sets."""
    n = len(H.elts)
    down = [set() for _ in range(n)]
    for a in range(n):
        stack = list(H.cov[a])
        while stack:
            x = stack.pop()
            if x not in down[a]:
                down[a].add(x)
                stack += list(H.cov[x])
    return [frozenset(d) for d in down]


CASES = [("D", 4, 2), ("D", 5, 2), ("D", 5, 3), ("A", 4, 2), ("A", 5, 3),
         ("D", 6, 3), ("E", 6, 2)]


def run(typ, n, k, cap=200, maxsize=13):
    rs = RootSystem(typ, n)
    lam = rs.fundamental(k)
    orb = orbit_with_words(rs, lam)
    if len(orb) > cap:
        return None
    good = bad = skipped = 0
    ex = []
    for mu, word in orb.items():
        if not word:
            continue
        H = build_heap(rs, word, lam)
        if not simple_spectrum(H) or len(H.elts) > maxsize:
            skipped += 1
            continue
        down = down_closure(H)
        R = rpp_series(down)
        P = hook_product(hooks(rs, mu))
        if R == P:
            good += 1
        else:
            bad += 1
            if len(ex) < 3:
                ex.append((mu, len(H.elts), hooks(rs, mu), R[:8], P[:8]))
    return good, bad, skipped, ex


if __name__ == "__main__":
    t0 = time.time()
    G = B = S = 0
    for typ, n, k in CASES:
        r = run(typ, n, k)
        if not r:
            continue
        g, b, s, ex = r
        G += g; B += b; S += s
        print(f"  {typ}{n} om{k}: {g}/{g+b} simple-spectrum elements match "
              f"({s} skipped)")
        for mu, m, H, R, P in ex:
            print(f"      mu={mu} |P|={m} hooks={H}")
            print(f"        RPP  {R}")
            print(f"        prod {P}")
    print(f"\n  hook product = RPP generating function: {G}/{G+B}")
    print(f"({time.time()-t0:.0f}s)")
    sys.exit(0 if B == 0 else 1)
