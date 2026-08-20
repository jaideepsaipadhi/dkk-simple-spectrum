"""
Direct computation of ch_q C[N_V] for G a product of copies of SL2, and a
DIRECT test of Theorem 1.2 -- not routed through Proposition 2.1.

Why this is the test that matters.  Proposition 2.1 derives Theorem 1.2 from
identity (2) only for COFREE V: it uses Lemma 2.3, C[N_V] (x) C[V]^G = C[V],
which needs cofreeness.  So a Hesselink-type V that is not cofree could satisfy
(2) and still violate Theorem 1.2, and a test based on (2) would never see it.
Here we instead build

        C[N_V] = C[V] / (C[V]^G_+ C[V])

degree by degree by linear algebra, decompose each graded piece into
irreducible characters, and compare against sum_lambda M_q(lambda) chi_lambda.

Setup for G = SL2^k.  V = sum_j V_{d_j} with d_j = (d_{j1},...,d_{jk}) and
V_d = tensor_i Sym^{d_i}(C^2).  A basis vector of the j-th summand is indexed
by (j, m_1..m_k) with 0 <= m_i <= d_{ji}, of weight (d_{j1}-2m_1, ...).  The
raising operator of the i-th sl2 acts by  e.(x^{d-m} y^m) = m x^{d-m+1} y^{m-1}.

On V* a coordinate xi_b satisfies X.xi_b = -sum_c M_{bc} xi_c where M is the
matrix of X on V; the action on C[V] is by derivations.  Invariants in degree
e are the weight-zero polynomials killed by every raising operator.

Linear algebra is done modulo a large prime; ranks are correct for all but
finitely many primes and the prime is reported.
"""
import sys, itertools
from collections import defaultdict

P = 1000003


# ------------------------------------------------------------ linear algebra
def rref(rows, ncols):
    mat = [r[:] for r in rows]
    piv, r = [], 0
    for c in range(ncols):
        p = next((i for i in range(r, len(mat)) if mat[i][c] % P), None)
        if p is None:
            continue
        mat[r], mat[p] = mat[p], mat[r]
        inv = pow(mat[r][c], P - 2, P)
        mat[r] = [(x * inv) % P for x in mat[r]]
        for i in range(len(mat)):
            if i != r and mat[i][c] % P:
                f = mat[i][c]
                mat[i] = [(mat[i][j] - f * mat[r][j]) % P for j in range(ncols)]
        piv.append(c)
        r += 1
        if r == len(mat):
            break
    return mat[:r], piv


def kernel(rows, ncols):
    R, piv = rref(rows, ncols)
    free = [c for c in range(ncols) if c not in piv]
    ker = []
    for f in free:
        v = [0] * ncols
        v[f] = 1
        for r, c in zip(R, piv):
            v[c] = (-r[f]) % P
        ker.append(v)
    return ker


# ------------------------------------------------------------------ the rep
class Rep:
    """V = sum_j V_{d_j} for G = SL2^k."""

    def __init__(self, summands):
        self.summands = [tuple(d) for d in summands]
        self.k = len(self.summands[0])
        self.basis = []                       # (j, m-tuple)
        for j, d in enumerate(self.summands):
            for m in itertools.product(*[range(x + 1) for x in d]):
                self.basis.append((j, m))
        self.index = {b: i for i, b in enumerate(self.basis)}
        self.dim = len(self.basis)

    def weight(self, b):
        j, m = b
        d = self.summands[j]
        return tuple(d[i] - 2 * m[i] for i in range(self.k))

    def raise_matrix(self, i):
        """M with e_i . v_c = sum_b M[b][c] v_b."""
        M = defaultdict(dict)
        for c, b in enumerate(self.basis):
            j, m = b
            if m[i] == 0:
                continue
            m2 = list(m); m2[i] -= 1
            tgt = self.index[(j, tuple(m2))]
            M[tgt][c] = m[i] % P
        return M


# ------------------------------------------------- monomials and the grading
def monomials(nvars, deg):
    if nvars == 0:
        yield ()
        return
    if nvars == 1:
        yield (deg,)
        return
    for a in range(deg + 1):
        for rest in monomials(nvars - 1, deg - a):
            yield (a,) + rest


def mono_weight(rep, mono):
    w = [0] * rep.k
    for b, e in enumerate(mono):
        if e:
            wb = rep.weight(rep.basis[b])
            for i in range(rep.k):
                w[i] -= e * wb[i]          # coordinates have weight -wt(v)
    return tuple(w)


def graded_monomials(rep, deg):
    out = defaultdict(list)
    for mono in monomials(rep.dim, deg):
        out[mono_weight(rep, mono)].append(mono)
    return out


# ------------------------------------------------------------- raising action
def apply_raise(rep, Mi, mono, coeff=1):
    """Derivation action of e_i on a monomial, as {monomial: coeff}."""
    out = defaultdict(int)
    for b, e in enumerate(mono):
        if e == 0:
            continue
        # e_i . xi_b = - sum_c M[b][c] xi_c   (b is the row index)
        row = Mi.get(b, {})
        for c, val in row.items():
            m2 = list(mono)
            m2[b] -= 1
            m2[c] += 1
            out[tuple(m2)] = (out[tuple(m2)] - coeff * e * val) % P
    return {k: v % P for k, v in out.items() if v % P}


def default_alphas(k):
    """omega-coordinates of the simple roots for SL2^k: alpha_i = 2 omega_i."""
    return [tuple(2 if t == i else 0 for t in range(k)) for i in range(k)]


def invariants(rep, deg, alphas=None, cache={}):
    """Basis of C[V]^G_deg as vectors over the weight-0 monomials.

    `alphas` gives the omega-coordinates of the simple roots, i.e. the weight
    shift produced by each raising operator.  For SL2^k this is 2*omega_i, but
    in general it is the i-th column of the Cartan matrix, so it must be passed
    in for any other group."""
    if alphas is None:
        alphas = default_alphas(rep.k)
    key = (tuple(rep.summands), deg, tuple(alphas))
    if key in cache:
        return cache[key]
    gm = graded_monomials(rep, deg)
    zero = tuple([0] * rep.k)
    base = gm.get(zero, [])
    if not base:
        cache[key] = ([], [])
        return cache[key]
    pos = {m: i for i, m in enumerate(base)}
    rows = []
    for i in range(rep.k):
        Mi = rep.raise_matrix(i)
        tgtw = tuple(alphas[i])
        tgt = {m: j for j, m in enumerate(gm.get(tgtw, []))}
        block = [[0] * len(base) for _ in range(len(tgt))]
        for m, col in pos.items():
            for m2, v in apply_raise(rep, Mi, m).items():
                if m2 in tgt:
                    block[tgt[m2]][col] = (block[tgt[m2]][col] + v) % P
        rows += block
    ker = kernel(rows, len(base)) if rows else \
        [[1 if x == y else 0 for x in range(len(base))] for y in range(len(base))]
    cache[key] = (ker, base)
    return cache[key]


# ----------------------------------------------------------- the nullcone
def nullcone_character(rep, maxdeg, verbose=False, alphas=None):
    """{degree: {weight: dim}} for C[N_V] = C[V] / (C[V]^G_+ C[V])."""
    out = {}
    invs = {}
    for e in range(1, maxdeg + 1):
        invs[e] = invariants(rep, e, alphas)
    for d in range(maxdeg + 1):
        gm = graded_monomials(rep, d)
        idealdim = defaultdict(int)
        for w, monos in gm.items():
            pos = {m: i for i, m in enumerate(monos)}
            rows = []
            for e in range(1, d + 1):
                ker, base = invs[e]
                if not ker:
                    continue
                lower = graded_monomials(rep, d - e)
                for h_w, hs in lower.items():
                    # weight of g*h is weight(g) + weight(h) = 0 + h_w
                    if h_w != w:
                        continue
                    for g in ker:
                        for h in hs:
                            row = [0] * len(monos)
                            for idx, cval in enumerate(g):
                                if cval % P == 0:
                                    continue
                                prod = tuple(a + b for a, b in zip(base[idx], h))
                                row[pos[prod]] = (row[pos[prod]] + cval) % P
                            if any(row):
                                rows.append(row)
            R, _ = rref(rows, len(monos)) if rows else ([], [])
            idealdim[w] = len(R)
        out[d] = {w: len(gm[w]) - idealdim[w] for w in gm
                  if len(gm[w]) - idealdim[w] > 0}
        if verbose:
            print(f"      deg {d}: dim C[V]_d = {sum(len(v) for v in gm.values())},"
                  f"  dim C[N_V]_d = {sum(out[d].values())}")
    return out


# --------------------------------------------- decompose a weight multiset
def decompose(wtmult, k):
    """Weight multiplicities -> irreducible multiplicities for SL2^k."""
    rem = dict(wtmult)
    out = defaultdict(int)
    while True:
        rem = {w: m for w, m in rem.items() if m}
        if not rem:
            break
        hw = max(rem, key=lambda w: (sum(w), w))
        if any(x < 0 for x in hw):
            raise ValueError(f"negative highest weight {hw}: not a character")
        mult = rem[hw]
        out[hw] += mult
        for m in itertools.product(*[range(0, x + 1) for x in hw]):
            w = tuple(hw[i] - 2 * m[i] for i in range(k))
            rem[w] = rem.get(w, 0) - mult
    return dict(out)
