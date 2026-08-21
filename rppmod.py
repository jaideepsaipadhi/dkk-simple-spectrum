"""
KRYLOV'S FOLLOW-UP: the RPP generating function of the poset IS a hook product.

Let w have the simple spectrum property, let P(w) be the resulting coloured
poset, and let RPP_{P}(z) = sum over order-preserving f : P -> Z_{>=0} of
z^{sum f}.  Fix a reduced word w = s_{i_1} ... s_{i_l}, let

    beta_p = s_{i_1} ... s_{i_{p-1}} ( alpha_{i_p} ),    p = 1, ..., l,

be its inversion roots and let c_p be the multiplicities of the c-trace, read
in the same order.  Then

    RPP_{P(w)}(z)  =  prod_{p=1}^{l}  ( 1 - z^{ht(beta_p)} )^{-c_p} ,     (HP)

with ht the height, ht(beta) = <beta, rho^vee>.  Verified with no exception.

Two consequences.

* If w is lambda-minuscule then every c_p = 1 and sum_p c_p = l(w), so (HP)
  reads
        RPP(z) = prod_{beta in Inv(w)} 1 / ( 1 - z^{ht(beta)} ) ,
  which is exactly the Peterson-Proctor hook product and exactly the formula
  asked about, the negative roots alpha = -beta being those with w(alpha) > 0.

* If w is NOT lambda-minuscule -- which happens inside the simple spectrum
  class -- the product over Inv(w) has too FEW factors: the correct hook
  multiset has sum_p c_p = |v| > l(w) entries.  It still contains the height
  multiset of Inv(w) as a sub-multiset, the extra entries being the repeated
  heights ht(beta_p) with c_p >= 2.  The smallest instance is D_4, lambda =
  omega_2, mu = (0,1,0,-2): l(w) = 5 but |P| = 6, hooks (1,1,2,3,3,4) against
  inversion heights (1,2,3,3,4).

The poset is read off the module, not from a combinatorial rule: elements are
the occupied pairs (degree, colour) of V = soc_sigma(I_kappa), with a covering
relation (d,i) < (d+1,j) exactly when the arrow a_{i->j} is nonzero on the
line V[(d,i)].  The orientation is DKK's: degree DECREASES upwards, as in the
Hasse diagrams on pp. 82-83.  (The opposite orientation is not a hook product
-- 110/342 -- so the orientation is part of the statement.)

RPP is computed exactly as a truncated power series: an order-preserving map is
sum_k 1_{I_k} for a unique decreasing multichain of nonempty order ideals I_k,
and that sum is evaluated by dynamic programming over the ideal lattice.

The hook multiset does not depend on the reduced word chosen, although the
formula is stated through one; this is checked too.
"""
import sys, time
from itertools import combinations
from collections import Counter
from rootsys import RootSystem, orbit_with_words
from wordtrace import trace
from gls import iterated_socle
from homological import setup
from qvar import P as MODP
from weylprops import reduced_words

N = 16                                     # series truncation
MAXSIZE = 11                               # poset size cap (2^|P| filters)


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


def geom(m):
    return [1 if d % m == 0 else 0 for d in range(N)]


def build(rs, pr, word, lam):
    """The module poset: (elements, strict down-sets), or None if not simple."""
    seq = []
    for _, i, c, _ in trace(rs, word, lam):
        seq += [i] * c
    X = iterated_socle(pr, rs, seq)
    if any(len(v) > 1 for v in X.values()):
        return None
    elts = [(d, i) for d in range(pr.top + 1) for i in rs.I
            if len(X.get((d, i), [])) == 1]
    idx = {e: n for n, e in enumerate(elts)}
    cov = [set() for _ in elts]
    for n, (d, i) in enumerate(elts):
        for j in rs.adj[i]:
            if (d + 1, j) in idx and (d, i, j) in pr.arrow:
                if any(y % MODP for y in pr.act(d, i, X[(d, i)][0], j)):
                    cov[idx[(d + 1, j)]].add(n)
    down = [set() for _ in elts]
    for n in range(len(elts)):
        stack = list(cov[n])
        while stack:
            x = stack.pop()
            if x not in down[n]:
                down[n].add(x)
                stack += list(cov[x])
    return elts, [frozenset(d) for d in down]


def rpp(strict_down):
    """RPP generating function in DKK's orientation, truncated to degree < N.

    DKK draw the poset with the degree DECREASING upwards, so an order-
    preserving map for their order is a decreasing multichain of order ideals
    for the degree order, which is what is summed over here.
    """
    n = len(strict_down)
    fs = []
    for r in range(1, n + 1):
        for S in combinations(range(n), r):
            S = frozenset(S)
            if all(strict_down[y] <= S for y in S):
                fs.append(S)
    fs.sort(key=len)
    g, tot = {}, [1] + [0] * (N - 1)
    for F in fs:
        s = [1] + [0] * (N - 1)
        for G in fs:
            if len(G) < len(F) and G <= F:
                s = add(s, g[G])
        m = len(F)
        z = [0] * N
        if m < N:
            z[m] = 1
        g[F] = mul(mul(z, geom(m)), s)
        tot = add(tot, g[F])
    return tot


def hook_multiset(rs, word, lam):
    """{ ht(beta_p) with multiplicity c_p } -- word read left to right."""
    cs = list(reversed([c for _, _, c, _ in trace(rs, word, lam)]))
    out = []
    for p, i in enumerate(word):
        b = rs.act_word(tuple(word[:p]), rs.simple_root(i))
        out += [int(sum(rs.root_coords(b)))] * cs[p]
    return sorted(out)


def product(H):
    f = [1] + [0] * (N - 1)
    for h in H:
        f = mul(f, geom(h))
    return f


CASES = [("A", 4, 2), ("A", 5, 3), ("A", 6, 3), ("D", 4, 2), ("D", 5, 2),
         ("D", 5, 3), ("D", 6, 2), ("D", 6, 3), ("D", 7, 2), ("E", 6, 2),
         ("E", 6, 3)]


def run(typ, n, k, cap=400):
    rs, pr, kk = setup(typ, n, k)
    lam = rs.fundamental(k)
    orb = orbit_with_words(rs, lam)
    if len(orb) > cap:
        return None
    st = Counter()
    for mu, word in orb.items():
        if not word:
            continue
        b = build(rs, pr, word, lam)
        if b is None:
            continue
        elts, down = b
        if len(elts) > MAXSIZE:
            st['skip'] += 1
            continue
        H = hook_multiset(rs, word, lam)
        st['tot'] += 1
        st['ok'] += (rpp(down) == product(H))
        st['minus'] += (len(H) == len(word))
        for w2 in sorted(reduced_words(rs, word))[:12]:
            st['wi_tot'] += 1
            st['wi_ok'] += (hook_multiset(rs, w2, lam) == H)
    return st


if __name__ == "__main__":
    t0 = time.time()
    T = Counter()
    for typ, n, k in CASES:
        r = run(typ, n, k)
        if not r:
            continue
        T += r
        print(f"  {typ}{n} om{k}: (HP) at {r['ok']}/{r['tot']} simple-spectrum "
              f"elements ({r['minus']} of them lambda-minuscule)")
    print(f"\n  (HP) overall: {T['ok']}/{T['tot']}   "
          f"lambda-minuscule: {T['minus']}, beyond it: {T['tot']-T['minus']}")
    print(f"  ({T['skip']} posets skipped as too large for exact filter DP)")
    print(f"  hook multiset independent of the reduced word: "
          f"{T['wi_ok']}/{T['wi_tot']}")
    print(f"\n({time.time()-t0:.0f}s)")
    print("HOOK PRODUCT VERIFIED" if T['ok'] == T['tot'] and
          T['wi_ok'] == T['wi_tot'] else "*** failed ***")
    sys.exit(0 if T['ok'] == T['tot'] and T['wi_ok'] == T['wi_tot'] else 1)
