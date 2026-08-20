"""Extend the submodule search to sums of two cyclic submodules."""
from itertools import product, combinations
from rootsys import RootSystem
from wordtrace import trace, v_from_trace
from qvar import Proj, rref, P
from qsearch import pieces, gen_submodule, quot_dims


def all_cyclic(pr, pc):
    out = []
    for (d0, i0), cols in sorted(pc.items()):
        n = len(cols)
        for coeffs in product(range(P), repeat=n):
            if all(a == 0 for a in coeffs): continue
            if next(a for a in coeffs if a) != 1: continue
            out.append(gen_submodule(pr, pc, d0, i0, coeffs))
    return out


def add(pr, A, B):
    keys = set(A) | set(B)
    return {k: rref(A.get(k, []) + B.get(k, []), len(pr.basis[k[0]]))[0]
            for k in keys}


def totals(rs, pr, pc, N):
    q = quot_dims(pr, pc, N)
    tot = {}
    for (d, i), c in q.items(): tot[i] = tot.get(i, 0) + c
    rows = {}
    for (d, i), c in q.items():
        if c: rows.setdefault(d, {})[i] = c
    return tot, rows


def search2(rs, k, v, dmax=8):
    pr = Proj(rs, k, dmax=dmax); pc = pieces(pr)
    target = dict(zip(rs.I, v))
    cyc = all_cyclic(pr, pc)
    found = {}
    def note(N):
        tot, rows = totals(rs, pr, pc, N)
        if all(tot.get(i, 0) == target.get(i, 0) for i in rs.I):
            sig = tuple(sorted((d, tuple(sorted(r.items()))) for d, r in rows.items()))
            found[sig] = rows
    for N in cyc: note(N)
    if not found:
        for A, B in combinations(cyc, 2):
            note(add(pr, A, B))
    return found


if __name__ == "__main__":
    for label, (typ, n, k), word in [
        ("Ex 5.10", ("D",4,2), (1,2,3,4,2)),
        ("Ex 5.11", ("D",5,2), (2,1,3,4,5,3,2)),
        ("Ex 5.12", ("D",5,2), (2,1,2,3,4,5,3,2)),
    ]:
        rs = RootSystem(typ, n)
        v = v_from_trace(rs, trace(rs, word, rs.fundamental(k)))
        f = search2(rs, k, v, dmax=8)
        print(f"{label}  v={v}:  {len(f)} distinct graded-dimension solution(s)")
        for sig, rows in f.items():
            for d in sorted(rows): print(f"      deg {d}: {rows[d]}")
            print(f"      simple spectrum: {all(c<=1 for r in rows.values() for c in r.values())}")
        print()
