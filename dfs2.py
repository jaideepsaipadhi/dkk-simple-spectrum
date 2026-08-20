"""
Optimised complete search.

Same top-down DFS over graded subspaces, plus:
  * trim dmax to the actual top nonzero degree of P
  * prune on remaining capacity: if acc + (all remaining pieces) < target, cut
  * prune per vertex as soon as acc exceeds target
  * optional early exit at the first profile (uniqueness verified on 192 cases)
  * cache subspace lists by dimension
"""
from itertools import product
from functools import lru_cache
from qvar import Proj, rref, P
from dfs import subspaces, Ambient

_SUB = {}
def subs(m):
    if m not in _SUB: _SUB[m] = subspaces(m)
    return _SUB[m]


def solve(rs, k, v, dmax=14, first_only=True):
    A = Ambient(rs, k, dmax)
    top = max((d for (d, i) in A.dim if A.dim[(d, i)]), default=0)
    target = dict(zip(rs.I, v))
    # capacity[d][i] = total dim of P[d'][i] for d' <= d
    cap = {}
    for i in rs.I:
        run = 0
        for d in range(top + 1):
            run += A.dim.get((d, i), 0)
            cap[(d, i)] = run
    results = {}

    def rec(d, N, acc):
        if results and first_only: return
        if d < 0:
            if all(acc.get(i, 0) == target.get(i, 0) for i in rs.I):
                rows = {}
                for (dd, ii), cols in A.cols.items():
                    q = len(cols) - len(N.get((dd, ii), []))
                    if q: rows.setdefault(dd, {})[ii] = q
                sig = tuple(sorted((x, tuple(sorted(y.items()))) for x, y in rows.items()))
                results[sig] = rows
            return
        # capacity prune
        for i in rs.I:
            if acc.get(i, 0) + cap.get((d, i), 0) < target.get(i, 0): return
        verts = [i for i in rs.I if A.dim.get((d, i), 0)]
        if not verts:
            rec(d - 1, N, acc); return
        pres = {i: A.preimage(d, i, N) for i in verts}
        # order vertices by fewest options first
        verts.sort(key=lambda i: len(subs(len(pres[i]))))
        def walk(t, N2, acc2):
            if results and first_only: return
            if t == len(verts):
                rec(d - 1, N2, acc2); return
            i = verts[t]; basis = pres[i]; m = A.dim[(d, i)]
            for sub in subs(len(basis)):
                q = m - len(sub)
                a = acc2.get(i, 0) + q
                if a > target.get(i, 0): continue
                if a + cap.get((d - 1, i), 0) < target.get(i, 0): continue
                rows = [[sum(co[t2] * basis[t2][c] for t2 in range(len(basis))) % P
                         for c in range(m)] for co in sub]
                R, _ = rref(rows, m)
                N3 = dict(N2); N3[(d, i)] = R
                acc3 = dict(acc2); acc3[i] = a
                walk(t + 1, N3, acc3)
        walk(0, N, acc)

    rec(top, {}, {})
    return results


if __name__ == "__main__":
    import time
    from rootsys import RootSystem
    from wordtrace import trace, v_from_trace
    for label, (typ, n, k), word in [
        ("Ex 5.10", ("D",4,2), (1,2,3,4,2)),
        ("Ex 5.11", ("D",5,2), (2,1,3,4,5,3,2)),
        ("Ex 5.12", ("D",5,2), (2,1,2,3,4,5,3,2)),
    ]:
        rs = RootSystem(typ, n)
        v = v_from_trace(rs, trace(rs, word, rs.fundamental(k)))
        t = time.time(); r = solve(rs, k, v, dmax=8)
        rows = list(r.values())[0]
        ok = all(c <= 1 for x in rows.values() for c in x.values())
        print(f"{label}: {dict(sorted(rows.items()))}  simple={ok}  ({time.time()-t:.2f}s)")
