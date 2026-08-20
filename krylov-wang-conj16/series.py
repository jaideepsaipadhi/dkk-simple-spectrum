"""
Hilbert series needed to test identity (2) of Krylov-Wang:

                            dim_q C[V^{rho^vee}]^h
                     X_V = ------------------------                        (2)
                                dim_q C[V]^G

Both sides are computed independently of X_V, so (2) is a real test.

  * dim_q C[V]^G is the multiplicity of the trivial character in
    ch_q C[V] = prod_{alpha in S} 1/(1 - q t^{-alpha}).  We extract it by
    multiplying by the Weyl denominator A(rho) = sum_w (-1)^l(w) t^{w rho}
    and reading the coefficient of t^{rho}: for any lambda,
        mult of chi_lambda in degree d  =  [t^{lambda+rho}] ( A(rho) * ch_q C[V] ).

  * dim_q C[V^{rho^vee}]^h counts monomials in the height-zero weight spaces
    whose total weight is 0, i.e. functions m : S0 -> Z_{>=0} with
    sum m(alpha) alpha = 0, graded by |m|.  A direct convolution.

Everything is truncated at order q^N.
"""
from fractions import Fraction as F
from collections import defaultdict
from rootdata import Root
from hess import rep_weights, split


def zero(N):
    return [0] * (N + 1)


def mulseries(a, b, N):
    """a, b: dict weight -> list of q-coefficients (truncated at N)."""
    out = defaultdict(lambda: zero(N))
    for wa, pa in a.items():
        for wb, pb in b.items():
            w = tuple(x + y for x, y in zip(wa, wb))
            row = out[w]
            for i, ca in enumerate(pa):
                if ca == 0:
                    continue
                for j, cb in enumerate(pb):
                    if cb and i + j <= N:
                        row[i + j] += ca * cb
    return {k: v for k, v in out.items() if any(v)}


def ch_poly_ring(rt, S, N):
    """prod_{alpha in S} 1/(1 - q t^{-alpha}), truncated at q^N."""
    cur = {tuple([F(0)] * rt.n): [1] + [0] * N}
    for a, mult in S.items():
        for _ in range(mult):
            fac = {}
            for k in range(N + 1):
                w = tuple(-k * x for x in a)
                row = fac.setdefault(w, zero(N))
                row[k] += 1
            cur = mulseries(cur, fac, N)
    return cur


def weyl_denominator(rt):
    out = {}
    for word, ln in rt.weyl():
        w = rt.apply_word(word, rt.rho)
        out[w] = out.get(w, 0) + (-1) ** ln
    return {k: [v] for k, v in out.items() if v}


def multiplicities(rt, S, N, lam=None):
    """q-series of the multiplicity of chi_lam (default: trivial) in C[V]."""
    if lam is None:
        lam = tuple([F(0)] * rt.n)
    ch = ch_poly_ring(rt, S, N)
    A = {k: v + [0] * N for k, v in weyl_denominator(rt).items()}
    prod = mulseries(A, ch, N)
    target = tuple(lam[k] + rt.rho[k] for k in range(rt.n))
    return prod.get(target, zero(N))


def invariants_series(rt, S, N):
    return multiplicities(rt, S, N)


def zero_weight_subring(rt, S0, N):
    """dim_q C[V^{rho^vee}]^h: functions m : S0 -> Z_{>=0}, sum m(a) a = 0."""
    cur = {tuple([F(0)] * rt.n): [1] + [0] * N}
    for a in S0:
        fac = {}
        for k in range(N + 1):
            w = tuple(k * x for x in a)
            row = fac.setdefault(w, zero(N))
            row[k] += 1
        cur = mulseries(cur, fac, N)
    z = tuple([F(0)] * rt.n)
    return cur.get(z, zero(N))


def polymul(a, b, N):
    out = zero(N)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                if y and i + j <= N:
                    out[i + j] += x * y
    return out


def poly_from_dict(p, N):
    out = zero(N)
    for e, c in p.items():
        if e <= N:
            out[e] += c
    return out


def test_identity2(rt, summands, N=10):
    """Check X_V * dim_q C[V]^G == dim_q C[V^{rho^vee}]^h up to q^N.

    Returns (holds, lhs, rhs, XV_poly).  Only meaningful when V is
    Hesselink-type, so that X_V is a scalar polynomial."""
    from hess import is_hesselink
    ok, X = is_hesselink(rt, summands)
    triv = tuple([F(0)] * rt.n)
    xp = poly_from_dict(X.get(triv, {}), N)
    S = rep_weights(rt, summands)
    Sp, S0, Sm = split(rt, S)
    inv = invariants_series(rt, S, N)
    zws = zero_weight_subring(rt, S0, N)
    lhs = polymul(xp, inv, N)
    return (lhs == zws), lhs, zws, xp, ok


if __name__ == "__main__":
    print("Validation on cases where Theorem 1.2 is proved\n")
    N = 8
    cases = [("A1 adjoint", Root("A", 1), [(2,)]),
             ("A2 adjoint", Root("A", 2), [(1, 1)]),
             ("B2 adjoint", Root("B", 2), [(0, 2)]),
             ("G2 adjoint", Root("G", 2), [(1, 0)])]
    for name, rt, sm in cases:
        holds, lhs, rhs, xp, hess = test_identity2(rt, sm, N)
        print(f"  {name}: Hesselink={hess}  identity (2) holds up to q^{N}: {holds}")
        print(f"      X_V              = {xp}")
        print(f"      dim_q C[V]^G     = {invariants_series(rt, rep_weights(rt, sm), N)}")
        print(f"      X_V * inv        = {lhs}")
        print(f"      dim_q C[V^rho]^h = {rhs}")


# --------------------------------------------------- coregularity (necessary
# condition for cofreeness)
def polynomial_ring_degrees(H, N):
    """If H is the Hilbert series of a polynomial ring, return the multiset of
    generator degrees; else return None.  Greedy: the lowest positive degree d
    with a positive coefficient must be a generator degree, with multiplicity
    that coefficient; divide out (1-q^d)^{a_d} and repeat.  Valid up to q^N."""
    cur = list(H[:N + 1])
    if not cur or cur[0] != 1:
        return None
    degs = []
    for d in range(1, N + 1):
        a = cur[d]
        if a < 0:
            return None
        for _ in range(a):
            degs.append(d)
            # multiply cur by (1 - q^d)
            new = list(cur)
            for i in range(N, d - 1, -1):
                new[i] -= cur[i - d]
            cur = new
        if len(degs) > 4 * N:
            return None
    return degs if all(c == 0 for c in cur[1:]) else None


def is_coregular(rt, summands, N=10):
    """Necessary condition for cofreeness: C[V]^G is a polynomial algebra."""
    S = rep_weights(rt, summands)
    H = invariants_series(rt, S, N)
    return polynomial_ring_degrees(H, N), H
