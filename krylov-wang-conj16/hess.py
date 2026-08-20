"""
The Hesselink-type test of Krylov-Wang (arXiv:2608.03314, Definition 1.5).

                 sum_w (-1)^l(w) t^{w rho} prod_{alpha in S+} (1 - q t^{-w alpha})
        X_V  =   ------------------------------------------------------------------
                              sum_w (-1)^l(w) t^{w rho}

and V is *Hesselink-type* when X_V is a Z[q]-multiple of the trivial character.

Expanding the product over subsets T of S+ and writing sigma_T for the sum of
the weights in T, the numerator is

        sum_T (-q)^{|T|} A(rho - sigma_T),      A(mu) := sum_w (-1)^l(w) t^{w mu}.

Since A(mu) = 0 when mu is W-singular and A(mu) = (-1)^{l(y)} A(nu) when
nu = y(mu) is dominant regular, the Weyl character formula gives

        X_V = sum_T (-q)^{|T|} * sign(T) * chi_{nu(T) - rho},

so X_V is computed exactly, as a finite Z[q]-combination of irreducible
characters, with no series truncation.  The subsets are summed by a
convolution over the weights of S+, so the cost is |S+| times the number of
distinct partial sums, not 2^{|S+|}.

V is Hesselink-type iff every coefficient except the one at the trivial
character vanishes.
"""
from fractions import Fraction as F
from collections import defaultdict
from rootdata import Root


# --------------------------------------------------------------- tiny Z[q]
def padd(a, b):
    out = dict(a)
    for k, v in b.items():
        out[k] = out.get(k, 0) + v
        if out[k] == 0:
            del out[k]
    return out


def pstr(p):
    if not p:
        return "0"
    terms = []
    for k in sorted(p):
        c = p[k]
        if k == 0:
            terms.append(f"{c:+d}")
        elif k == 1:
            terms.append(f"{c:+d}q" if abs(c) != 1 else ("+q" if c > 0 else "-q"))
        else:
            terms.append(f"{c:+d}q^{k}" if abs(c) != 1
                         else (f"+q^{k}" if c > 0 else f"-q^{k}"))
    return "".join(terms).lstrip("+")


# ----------------------------------------------------------------- the rep
def rep_weights(rt, summands):
    """summands: list of highest weights (tuples).  Returns the multiset S."""
    S = defaultdict(int)
    for hw in summands:
        for w, m in rt.weights(hw).items():
            S[w] += m
    return dict(S)


def split(rt, S):
    Sp, S0, Sm = [], [], []
    for w, m in S.items():
        h = rt.height(w)
        tgt = Sp if h > 0 else (S0 if h == 0 else Sm)
        tgt += [w] * m
    return Sp, S0, Sm


# --------------------------------------------------------------------- X_V
def XV(rt, summands, verbose=False):
    """Return {dominant weight: Z[q] coefficient} for X_V."""
    S = rep_weights(rt, summands)
    Sp, S0, Sm = split(rt, S)
    # convolution: dict sigma -> Z[q]
    cur = {tuple([F(0)] * rt.n): {0: 1}}
    for a in Sp:
        nxt = defaultdict(dict)
        for sig, poly in cur.items():
            nxt[sig] = padd(nxt[sig], poly)
            sig2 = tuple(sig[k] + a[k] for k in range(rt.n))
            nxt[sig2] = padd(nxt[sig2], {e + 1: -c for e, c in poly.items()})
        cur = {k: v for k, v in nxt.items() if v}
    out = defaultdict(dict)
    rho = rt.rho
    for sig, poly in cur.items():
        mu = tuple(rho[k] - sig[k] for k in range(rt.n))
        nu, sign = rt.dominant_conjugate(mu)
        if any(x == 0 for x in nu):        # singular: A(mu) = 0
            continue
        lam = tuple(nu[k] - rho[k] for k in range(rt.n))
        out[lam] = padd(out[lam], {e: sign * c for e, c in poly.items()})
    return {k: v for k, v in out.items() if v}


def is_hesselink(rt, summands):
    X = XV(rt, summands)
    triv = tuple([F(0)] * rt.n)
    return all(k == triv for k in X), X


def show(name, rt, summands):
    ok, X = is_hesselink(rt, summands)
    triv = tuple([F(0)] * rt.n)
    bits = []
    for lam in sorted(X, key=lambda z: (sum(z), z)):
        tag = "1" if lam == triv else "chi" + str(tuple(int(x) for x in lam))
        bits.append(f"({pstr(X[lam])})*{tag}")
    print(f"  {name}: {'HESSELINK-TYPE' if ok else 'not Hesselink'}")
    print(f"      X_V = {' + '.join(bits) if bits else '0'}")
    return ok


if __name__ == "__main__":
    print("Reproducing the paper's examples\n")

    print("Theorem 1.2 case 1: the adjoint representation should be Hesselink-type")
    for t, n, hw in [("A", 1, (2,)), ("A", 2, (1, 1)), ("B", 2, (0, 2)),
                     ("G", 2, (1, 0))]:
        rt = Root(t, n)
        show(f"{t}{n} adjoint", rt, [hw])

    print("\nExample 4.2: sl2, V = V_3 -- cofree but NOT Hesselink-type")
    show("A1 V3", Root("A", 1), [(3,)])

    print("\nExample 4.3: sl2 x sl2, V = V_1 box V_3 -- NOT Hesselink-type")
    print("           (paper: X_V = 1 + q chi_{V1 box V1} - q^2 chi_{V0 box V4})")
    show("A1xA1 V1|V3", Root([("A", 1), ("A", 1)]), [(1, 3)])
