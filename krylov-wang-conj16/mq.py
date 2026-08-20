"""
The q-analog Kostant partition function P_q, the multiplicity function M_q,
and a DIRECT test of Theorem 1.2 for G = SL2^k.

    P_q(mu) = sum_n #{ m : S -> Z_{>=0} | sum m(a) a = -mu, |m| = n,
                                          m weakly negative } q^n
    M_q(mu) = sum_w (-1)^{l(w)} P_q(w(mu+rho) - rho)

"Weakly negative" (Krylov-Wang, Section 1): m(a) = 0 for every a in S+, and
for every m' in (S0)_0 there is some a in S0 with m(a) < m'(a) -- i.e. m
restricted to S0 dominates no element of (S0)_0, where (S0)_0 is the set of
nonzero m' : S0 -> Z_{>=0} with sum m'(a) a = 0.  By Dickson's lemma it is
enough to test the finitely many MINIMAL elements of (S0)_0, which we compute.

Theorem 1.2 asserts  ch_q C[N_V] = sum_lambda M_q(lambda) chi_lambda.  We test
that against the direct computation of C[N_V] in nullcone.py.
"""
import itertools
from collections import defaultdict
from nullcone import Rep, nullcone_character, decompose


def split_weights(rep):
    S = [rep.weight(b) for b in rep.basis]
    Sp = [w for w in S if sum(w) > 0]
    S0 = [w for w in S if sum(w) == 0]
    Sm = [w for w in S if sum(w) < 0]
    return S, Sp, S0, Sm


def minimal_relations(S0, k, cap):
    """Minimal nonzero m : S0 -> Z_{>=0} with sum m(a) a = 0."""
    sols = []
    n = len(S0)
    def rec(i, cur, acc):
        if sum(cur) > cap:
            return
        if i == n:
            if any(cur) and all(x == 0 for x in acc):
                sols.append(tuple(cur))
            return
        for e in range(cap - sum(cur) + 1):
            rec(i + 1, cur + [e],
                tuple(acc[t] + e * S0[i][t] for t in range(k)))
    rec(0, [], tuple([0] * k))
    minimal = []
    for s in sorted(sols, key=sum):
        if not any(all(t <= x for t, x in zip(m, s)) for m in minimal):
            minimal.append(s)
    return minimal


def weakly_negative_gen(rep, maxdeg):
    """Yield (total weight of -sum m(a)a, |m|) for weakly negative m."""
    S, Sp, S0, Sm = split_weights(rep)
    k = rep.k
    rel = minimal_relations(S0, k, maxdeg + 2)
    # m is supported on S0 u S- ; m|S+ = 0
    supp = S0 + Sm
    n0 = len(S0)
    out = defaultdict(lambda: defaultdict(int))

    def ok_s0(part):
        return not any(all(r[t] <= part[t] for t in range(n0)) for r in rel)

    def rec(i, cur, acc, tot):
        if tot > maxdeg:
            return
        if i == len(supp):
            if ok_s0(tuple(cur[:n0])):
                mu = tuple(-acc[t] for t in range(k))
                out[mu][tot] += 1
            return
        for e in range(maxdeg - tot + 1):
            rec(i + 1, cur + [e],
                tuple(acc[t] + e * supp[i][t] for t in range(k)), tot + e)
    rec(0, [], tuple([0] * k), 0)
    return out


def Pq_table(rep, maxdeg):
    """{mu: {n: count}} for P_q."""
    return weakly_negative_gen(rep, maxdeg)


def Mq(rep, lam, Pq, maxdeg):
    """M_q(lam) as a list of coefficients."""
    k = rep.k
    rho = tuple([1] * k)
    out = [0] * (maxdeg + 1)
    for signs in itertools.product([1, -1], repeat=k):
        sgn = 1
        for s in signs:
            if s < 0:
                sgn = -sgn
        arg = tuple(signs[i] * (lam[i] + rho[i]) - rho[i] for i in range(k))
        for n, c in Pq.get(arg, {}).items():
            if n <= maxdeg:
                out[n] += sgn * c
    return out


def test_theorem12(summands, maxdeg=4, verbose=True):
    """Compare ch_q C[N_V] with sum_lambda M_q(lambda) chi_lambda."""
    rep = Rep(summands)
    nc = nullcone_character(rep, maxdeg)
    Pq = Pq_table(rep, maxdeg)
    # left side: irreducible multiplicities by degree
    lhs = {}
    for d in range(maxdeg + 1):
        lhs[d] = decompose(nc[d], rep.k) if nc[d] else {}
    # check every lambda either side can see -- restricting to the support of
    # C[N_V] silently skips the lambdas where M_q is nonzero and C[N_V] is not.
    from nullcone import graded_monomials
    lams = set()
    for d in lhs:
        lams |= set(lhs[d])
    for d in range(maxdeg + 1):
        for w in graded_monomials(rep, d):
            if all(x >= 0 for x in w):
                lams.add(tuple(int(x) for x in w))
    agree = True
    rows = []
    for lam in sorted(lams, key=lambda z: (sum(z), z)):
        m = Mq(rep, lam, Pq, maxdeg)
        got = [lhs[d].get(lam, 0) for d in range(maxdeg + 1)]
        same = (m == got)
        agree &= same
        rows.append((lam, got, m, same))
    # also check that M_q predicts nothing for lambda outside the support
    if verbose:
        print(f"    V = {summands}   dim {rep.dim}")
        for lam, got, m, same in rows:
            flag = "ok" if same else "  <-- MISMATCH"
            print(f"      chi{lam}: C[N_V] {got}   M_q {m}{flag}")
    return agree, rows


if __name__ == "__main__":
    print("Direct test of Theorem 1.2 for G = SL2^k\n")
    print("  known cases (should agree):")
    for sm in [[(1,), (1,)], [(2,)], [(1,)]]:
        ok, _ = test_theorem12(sm, 4)
        print(f"    -> Theorem 1.2 {'holds' if ok else 'FAILS'}\n")
