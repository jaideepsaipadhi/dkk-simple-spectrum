"""
Root systems, Weyl groups and weight multiplicities for the Krylov-Wang
computation (arXiv:2608.03314).

Conventions.  Weights are stored in the basis of FUNDAMENTAL WEIGHTS, as
tuples of Fractions: lambda = sum a_i omega_i, so <lambda, alpha_i^vee> = a_i
and the simple reflection is

        s_i(lambda) = lambda - a_i * alpha_i.

The Cartan matrix is A[i][j] = <alpha_j, alpha_i^vee>.  In the omega-basis
alpha_j has coordinates (A[0][j], ..., A[r-1][j]).

Height.  rho^vee is the sum of the fundamental coweights, so <alpha_i,rho^vee>
= 1 for every simple root; hence for lambda = sum_j c_j alpha_j we have
<lambda, rho^vee> = sum_j c_j.  Passing from omega- to alpha-coordinates is
c = A^{-1} a.  This is `height`.
"""
from fractions import Fraction as F
from itertools import product


# ----------------------------------------------------------------- Cartan data
def cartan(typ, n):
    A = [[0] * n for _ in range(n)]
    for i in range(n):
        A[i][i] = 2
    def link(i, j, aij, aji):
        A[i][j] = aij      # <alpha_j, alpha_i^vee>
        A[j][i] = aji
    if typ == "A":
        for i in range(n - 1):
            link(i, i + 1, -1, -1)
    elif typ == "B":            # alpha_n short
        for i in range(n - 1):
            link(i, i + 1, -1, -1)
        if n >= 2:
            A[n - 2][n - 1] = -1      # alpha_n is SHORT in type B
            A[n - 1][n - 2] = -2
    elif typ == "C":            # alpha_n long
        for i in range(n - 1):
            link(i, i + 1, -1, -1)
        if n >= 2:
            A[n - 2][n - 1] = -2
            A[n - 1][n - 2] = -1
    elif typ == "D":
        for i in range(n - 2):
            link(i, i + 1, -1, -1)
        link(n - 3, n - 1, -1, -1)
    elif typ == "E":
        edges = {6: [(0, 2), (2, 3), (1, 3), (3, 4), (4, 5)],
                 7: [(0, 2), (2, 3), (1, 3), (3, 4), (4, 5), (5, 6)],
                 8: [(0, 2), (2, 3), (1, 3), (3, 4), (4, 5), (5, 6), (6, 7)]}[n]
        for i, j in edges:
            link(i, j, -1, -1)
    elif typ == "F":
        link(0, 1, -1, -1); A[1][2] = -1; A[2][1] = -2; link(2, 3, -1, -1)
    elif typ == "G":
        A[0][1] = -1; A[1][0] = -3     # alpha_1 long, alpha_2 short
    else:
        raise ValueError(typ)
    return A


def sym_factors(typ, n):
    """d_i = (alpha_i,alpha_i)/2, normalised so that short roots have d = 1."""
    if typ in "ADE":
        return [1] * n
    if typ == "B":
        return [2] * (n - 1) + [1]
    if typ == "C":
        return [1] * (n - 1) + [2]
    if typ == "F":
        return [2, 2, 1, 1]
    if typ == "G":
        return [3, 1]
    raise ValueError(typ)


def matinv(A):
    n = len(A)
    M = [[F(A[i][j]) for j in range(n)] + [F(1 if i == j else 0) for j in range(n)]
         for i in range(n)]
    for c in range(n):
        p = next(r for r in range(c, n) if M[r][c] != 0)
        M[c], M[p] = M[p], M[c]
        pv = M[c][c]
        M[c] = [x / pv for x in M[c]]
        for r in range(n):
            if r != c and M[r][c] != 0:
                f = M[r][c]
                M[r] = [M[r][k] - f * M[c][k] for k in range(2 * n)]
    return [row[n:] for row in M]


def block(components):
    """Cartan matrix and symmetrisers of a product of simple factors."""
    blocks = [(cartan(t, n), sym_factors(t, n)) for t, n in components]
    N = sum(len(b[0]) for b in blocks)
    A = [[0] * N for _ in range(N)]
    d = []
    off = 0
    for B, db in blocks:
        m = len(B)
        for i in range(m):
            for j in range(m):
                A[off + i][off + j] = B[i][j]
        d += db
        off += m
    return A, d


class Root:
    def __init__(self, typ, n=None):
        if isinstance(typ, (list, tuple)) and not isinstance(typ, str):
            self.components = list(typ)
            self.A, self.d = block(self.components)
            self.typ = "x".join(f"{t}{k}" for t, k in self.components)
            self.n = len(self.A)
        else:
            self.components = [(typ, n)]
            self.typ, self.n = typ, n
            self.A = cartan(typ, n)
            self.d = sym_factors(typ, n)
        n = self.n
        self.Ainv = matinv(self.A)
        for i in range(n):
            for j in range(n):
                assert self.d[i] * self.A[i][j] == self.d[j] * self.A[j][i], \
                    f"Cartan matrix of {typ}{n} is not symmetrisable by d={self.d}"

        self.rho = tuple(F(1) for _ in range(n))
        self.alpha = [tuple(F(self.A[i][j]) for i in range(n)) for j in range(n)]
        self._weyl = None
        self._posroots = None

    # ---- basic operations in omega-coordinates
    def refl(self, i, lam):
        a = lam[i]
        return tuple(lam[k] - a * self.alpha[i][k] for k in range(self.n))

    def height(self, lam):
        """<lambda, rho^vee> = sum of alpha-coordinates."""
        return sum(sum(self.Ainv[j][k] * lam[k] for k in range(self.n))
                   for j in range(self.n))

    def ip(self, lam, mu):
        """W-invariant inner product, normalised by d_i."""
        cl = [sum(self.Ainv[j][k] * lam[k] for k in range(self.n)) for j in range(self.n)]
        cm = [sum(self.Ainv[j][k] * mu[k] for k in range(self.n)) for j in range(self.n)]
        tot = F(0)
        for i in range(self.n):
            for j in range(self.n):
                tot += cl[i] * cm[j] * self.d[i] * self.A[i][j]
        return tot

    def is_dominant(self, lam):
        return all(x >= 0 for x in lam)

    # ---- Weyl group as a list of (permutation-of-weights closure) with lengths
    def weyl(self):
        """Return list of (word, length) and a callable action; built by BFS on
        the orbit of a regular dominant weight, which faithfully indexes W."""
        if self._weyl is not None:
            return self._weyl
        reg = self.rho
        seen = {reg: ()}
        frontier = [reg]
        while frontier:
            nxt = []
            for x in frontier:
                for i in range(self.n):
                    y = self.refl(i, x)
                    if y not in seen:
                        seen[y] = (i,) + seen[x]
                        nxt.append(y)
            frontier = nxt
        self._weyl = [(w, len(w)) for w in seen.values()]
        return self._weyl

    def apply_word(self, word, lam):
        for i in word:
            lam = self.refl(i, lam)
        return lam

    def positive_roots(self):
        if self._posroots is not None:
            return self._posroots
        seen = set(self.alpha)
        frontier = list(self.alpha)
        while frontier:
            nxt = []
            for x in frontier:
                for i in range(self.n):
                    y = self.refl(i, x)
                    if self.height(y) > 0 and y not in seen:
                        seen.add(y); nxt.append(y)
            frontier = nxt
        self._posroots = sorted(seen, key=lambda r: (self.height(r), r))
        return self._posroots

    # ---- make a weight dominant, recording the length of the element used
    def dominant_conjugate(self, lam):
        """Return (mu, sign) with mu dominant and mu = w(lam), sign = (-1)^l(w);
        returns (None, 0) if lam is on a wall (i.e. some coordinate hits 0
        during the reduction, meaning lam is W-conjugate to a non-regular
        weight)."""
        sign = 1
        lam = tuple(lam)
        for _ in range(10000):
            i = next((k for k in range(self.n) if lam[k] < 0), None)
            if i is None:
                return lam, sign
            lam = self.refl(i, lam)
            sign = -sign
        raise RuntimeError("no convergence")

    # ---- Freudenthal weight multiplicities
    def weights(self, hw):
        """Multiset of weights of L(hw) as {weight: multiplicity}."""
        hw = tuple(F(x) for x in hw)
        pos = self.positive_roots()
        # generate the weight system: all dominant weights <= hw, then W-orbits
        doms = self._dominant_below(hw)
        mult = {}
        rho = self.rho
        n2 = self.ip(hw, hw) + 2 * self.ip(hw, rho)
        for mu in doms:                      # in decreasing height order
            if mu == hw:
                mult[mu] = 1
                continue
            s = F(0)
            for a in pos:
                k = 1
                while True:
                    nu = tuple(mu[i] + k * a[i] for i in range(self.n))
                    dom, _ = self.dominant_conjugate(nu)
                    if dom not in mult:
                        break
                    s += mult[dom] * self.ip(nu, a)
                    k += 1
            den = n2 - self.ip(mu, mu) - 2 * self.ip(mu, rho)
            mult[mu] = F(2) * s / den if den != 0 else F(0)
        out = {}
        for mu, m in mult.items():
            for nu in self._orbit(mu):
                out[nu] = int(m)
        return out

    def _dominant_below(self, hw):
        """Dominant weights of the form hw - (nonneg combination of simple
        roots), sorted by decreasing height."""
        out, seen = [], set()
        frontier = [hw]
        seen.add(hw)
        while frontier:
            nxt = []
            for x in frontier:
                for i in range(self.n):
                    y = tuple(x[k] - self.alpha[i][k] for k in range(self.n))
                    if y in seen:
                        continue
                    dom, _ = self.dominant_conjugate(y)
                    seen.add(y)
                    if dom == y:
                        nxt.append(y)
                    else:
                        nxt.append(y)     # keep exploring through non-dominant
            frontier = [z for z in nxt if self.height(z) >= -self.height(hw)]
            out += frontier
        doms = [hw] + [x for x in out if self.is_dominant(x)]
        uniq = []
        s = set()
        for x in doms:
            if x not in s:
                s.add(x); uniq.append(x)
        uniq.sort(key=lambda m: -self.height(m))
        return uniq

    def _orbit(self, lam):
        seen = {tuple(lam)}
        frontier = [tuple(lam)]
        while frontier:
            nxt = []
            for x in frontier:
                for i in range(self.n):
                    y = self.refl(i, x)
                    if y not in seen:
                        seen.add(y); nxt.append(y)
            frontier = nxt
        return seen
