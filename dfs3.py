"""Complete graded-submodule search on the walk-free projective (proj2)."""
from qvar import rref, P
from proj2 import Proj2
from dfs import subspaces

_SUB = {}
def subs(m):
    if m not in _SUB: _SUB[m] = subspaces(m)
    return _SUB[m]


def solve(rs, k, v, first_only=True, _cache={}):
    key = (rs.typ, rs.n, k)
    if key not in _cache: _cache[key] = Proj2(rs, k)
    pr = _cache[key]
    top = pr.top
    target = dict(zip(rs.I, v))
    cap = {}
    for i in rs.I:
        run = 0
        for d in range(top + 1):
            run += pr.dim.get((d, i), 0)
            cap[(d, i)] = run
    results = {}

    def preimage(d, i, N):
        m = pr.dim.get((d, i), 0)
        if m == 0: return []
        eqs = []
        for j in rs.I:
            nj = pr.dim.get((d + 1, j), 0)
            if nj == 0 or (d, i, j) not in pr.arrow: continue
            R, piv = rref([r[:] for r in N.get((d + 1, j), [])], nj)
            imgs = []
            for b in range(m):
                e = [1 if t == b else 0 for t in range(m)]
                img = pr.act(d, i, e, j)
                for r, c in zip(R, piv):
                    if img[c] % P:
                        g = img[c]
                        img = [(img[t] - g * r[t]) % P for t in range(nj)]
                imgs.append(img)
            for t in range(nj):
                row = [imgs[b][t] % P for b in range(m)]
                if any(x % P for x in row): eqs.append(row)
        if not eqs:
            return [[1 if t == b else 0 for t in range(m)] for b in range(m)]
        R, piv = rref(eqs, m)
        free = [c for c in range(m) if c not in piv]
        ker = []
        for f in free:
            x = [0] * m; x[f] = 1
            for r, c in zip(R, piv): x[c] = (-r[f]) % P
            ker.append(x)
        return ker

    def rec(d, N, acc):
        if results and first_only: return
        if d < 0:
            if all(acc.get(i, 0) == target.get(i, 0) for i in rs.I):
                rows = {}
                for (dd, ii), m in pr.dim.items():
                    q = m - len(N.get((dd, ii), []))
                    if q: rows.setdefault(dd, {})[ii] = q
                sig = tuple(sorted((x, tuple(sorted(y.items()))) for x, y in rows.items()))
                results[sig] = rows
            return
        for i in rs.I:
            if acc.get(i, 0) + cap.get((d, i), 0) < target.get(i, 0): return
        verts = [i for i in rs.I if pr.dim.get((d, i), 0)]
        if not verts:
            rec(d - 1, N, acc); return
        pres = {i: preimage(d, i, N) for i in verts}
        def walk(t, N2, acc2):
            if results and first_only: return
            if t == len(verts):
                rec(d - 1, N2, acc2); return
            i = verts[t]; basis = pres[i]; m = pr.dim[(d, i)]
            for sub in subs(len(basis)):
                q = m - len(sub)
                a = acc2.get(i, 0) + q
                if a > target.get(i, 0): continue
                if a + cap.get((d - 1, i), 0) < target.get(i, 0): continue
                rows = [[sum(co[u] * basis[u][c] for u in range(len(basis))) % P
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
    for label, (typ,n,k), word in [
        ("Ex 5.10", ("D",4,2), (1,2,3,4,2)),
        ("Ex 5.11", ("D",5,2), (2,1,3,4,5,3,2)),
        ("Ex 5.12", ("D",5,2), (2,1,2,3,4,5,3,2)),
    ]:
        rs = RootSystem(typ, n)
        v = v_from_trace(rs, trace(rs, word, rs.fundamental(k)))
        t=time.time(); r = solve(rs, k, v)
        rows = list(r.values())[0]
        print(f"{label}: {dict(sorted(rows.items()))} "
              f"simple={all(c<=1 for x in rows.values() for c in x.values())} ({time.time()-t:.3f}s)")
