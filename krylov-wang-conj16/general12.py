"""
The direct Theorem 1.2 test for any product of type-A factors.

Generalises mq.py (which was SL2^k only) by using
  * genrep.GenRep      -- explicit representations, hence C[V]^G and C[N_V],
  * rootdata.Root      -- the Weyl group with lengths, for M_q,
  * a general "peel the highest weight" decomposition.

The SL2^k search found Hesselink-type = trivial-type exactly, i.e. exactly the
class the paper can already prove, so any counterexample to Conjecture 1.6 must
sit at another semisimple group.  This module is what tests those.
"""
import itertools, sys, time
from fractions import Fraction as F
from collections import defaultdict
from rootdata import Root
from genrep import GenRep
from nullcone import (graded_monomials, invariants, nullcone_character, rref, P)


def int_weight(w):
    return tuple(int(x) for x in w)


def rep_weight_multiset(rt, summands):
    S = defaultdict(int)
    for hw in summands:
        for w, m in rt.weights(tuple(F(x) for x in hw)).items():
            S[int_weight(w)] += m
    return dict(S)


def split_weights(rt, S):
    Sp, S0, Sm = [], [], []
    for w, m in S.items():
        h = rt.height(tuple(F(x) for x in w))
        tgt = Sp if h > 0 else (S0 if h == 0 else Sm)
        tgt += [w] * m
    return Sp, S0, Sm


def minimal_relations(S0, k, cap):
    sols, n = [], len(S0)
    def rec(i, cur, acc):
        if sum(cur) > cap:
            return
        if i == n:
            if any(cur) and all(x == 0 for x in acc):
                sols.append(tuple(cur))
            return
        for e in range(cap - sum(cur) + 1):
            rec(i + 1, cur + [e], tuple(acc[t] + e * S0[i][t] for t in range(k)))
    rec(0, [], tuple([0] * k))
    minimal = []
    for s in sorted(sols, key=sum):
        if not any(all(t <= x for t, x in zip(m, s)) for m in minimal):
            minimal.append(s)
    return minimal


def Pq_table(rt, S, maxdeg):
    k = rt.n
    Sp, S0, Sm = split_weights(rt, S)
    rel = minimal_relations(S0, k, maxdeg + 2)
    supp, n0 = S0 + Sm, len(S0)
    out = defaultdict(lambda: defaultdict(int))

    def ok_s0(part):
        return not any(all(r[t] <= part[t] for t in range(n0)) for r in rel)

    def rec(i, cur, acc, tot):
        if tot > maxdeg:
            return
        if i == len(supp):
            if ok_s0(tuple(cur[:n0])):
                out[tuple(-acc[t] for t in range(k))][tot] += 1
            return
        for e in range(maxdeg - tot + 1):
            rec(i + 1, cur + [e],
                tuple(acc[t] + e * supp[i][t] for t in range(k)), tot + e)
    rec(0, [], tuple([0] * k), 0)
    return out


def Mq(rt, lam, Pq, maxdeg):
    out = [0] * (maxdeg + 1)
    rho = rt.rho
    lamrho = tuple(F(lam[i]) + rho[i] for i in range(rt.n))
    for word, ln in rt.weyl():
        arg = rt.apply_word(word, lamrho)
        arg = tuple(arg[i] - rho[i] for i in range(rt.n))
        key = int_weight(arg)
        for n, c in Pq.get(key, {}).items():
            if n <= maxdeg:
                out[n] += ((-1) ** ln) * c
    return out


def decompose(rt, wtmult):
    """Weight multiplicities -> irreducible multiplicities, general type."""
    rem = {k: v for k, v in wtmult.items() if v}
    out = defaultdict(int)
    guard = 0
    while rem:
        guard += 1
        if guard > 5000:
            raise RuntimeError("decomposition did not terminate")
        doms = [w for w in rem if all(x >= 0 for x in w)]
        if not doms:
            raise ValueError(f"no dominant weight left: {rem}")
        hw = max(doms, key=lambda w: (rt.height(tuple(F(x) for x in w)), w))
        mult = rem[hw]
        if mult < 0:
            raise ValueError(f"negative multiplicity at {hw}")
        out[hw] += mult
        for w, m in rt.weights(tuple(F(x) for x in hw)).items():
            key = int_weight(w)
            rem[key] = rem.get(key, 0) - mult * m
        rem = {k: v for k, v in rem.items() if v}
    return dict(out)


def test(ranks, summands, maxdeg=3, verbose=True):
    rt = Root([("A", n) for n in ranks])
    rep = GenRep(ranks, summands)
    S = rep_weight_multiset(rt, summands)
    assert sum(S.values()) == rep.dim, "weight count disagrees with the model"
    alphas = [tuple(int(rt.A[t][i]) for t in range(rt.n)) for i in range(rt.n)]
    nc = nullcone_character(rep, maxdeg, alphas=alphas)
    Pq = Pq_table(rt, S, maxdeg)
    lhs = {d: (decompose(rt, nc[d]) if nc[d] else {}) for d in range(maxdeg + 1)}
    # IMPORTANT: check every lambda that either side can see.  Restricting to
    # the support of C[N_V] would silently skip the lambdas where M_q predicts
    # a nonzero multiplicity and C[N_V] has none -- exactly the asymmetry a
    # counterexample is likely to show.
    lams = set()
    for d in lhs:
        lams |= set(lhs[d])
    for d in range(maxdeg + 1):
        for w in graded_monomials(rep, d):
            if all(x >= 0 for x in w):
                lams.add(tuple(int(x) for x in w))
    agree, rows = True, []
    for lam in sorted(lams, key=lambda z: (sum(z), z)):
        m = Mq(rt, lam, Pq, maxdeg)
        got = [lhs[d].get(lam, 0) for d in range(maxdeg + 1)]
        same = (m == got)
        agree &= same
        rows.append((lam, got, m, same))
    if verbose:
        print(f"    ranks {ranks}  V = {summands}  dim {rep.dim}")
        for lam, got, m, same in rows:
            print(f"      chi{lam}: C[N_V] {got}   M_q {m}"
                  f"{'' if same else '   <-- MISMATCH'}")
    return agree, rows


if __name__ == "__main__":
    print("Calibration of the general test\n")
    print("  should HOLD (Theorem 1.2 case 1: the adjoint representation):")
    ok, _ = test([2], [(1, 1)], 3)
    print(f"    -> {'holds' if ok else 'FAILS'}\n")
    print("  should HOLD (SL2, V_1 + V_1):")
    ok, _ = test([1], [(1,), (1,)], 3)
    print(f"    -> {'holds' if ok else 'FAILS'}\n")
    print("  should FAIL (SL2, V_3 -- Section 2.1 gives it a modified P):")
    ok, _ = test([1], [(3,)], 3)
    print(f"    -> {'holds' if ok else 'FAILS'}")
