"""
Complete submodule search: top-down DFS over graded subspaces.

A graded submodule N of P = Pi e_k is a choice of subspace N[d][i] <= P[d][i]
with  act(N[d][i]) <= N[d+1][j]  for every arrow.  Sweeping degrees from the
TOP down, the constraint at degree d depends only on N[d+1], which is already
fixed, so each piece is an independent choice of subspace of a computable
preimage.  This finds EVERY graded submodule -- no generator-count limit.
"""
from itertools import product
from qvar import Proj, rref, P, inv


def subspaces(m):
    """All subspaces of F_P^m, as lists of rref rows."""
    if m == 0:
        return [[]]
    out = []
    for r in range(m + 1):
        for piv in _combos(range(m), r):
            free = [c for c in range(m) if c not in piv]
            slots = [(a, c) for a in range(r) for c in free if c > piv[a]]
            for vals in product(range(P), repeat=len(slots)):
                rows = [[0] * m for _ in range(r)]
                for a in range(r):
                    rows[a][piv[a]] = 1
                for (a, c), val in zip(slots, vals):
                    rows[a][c] = val
                out.append(rows)
    return out


def _combos(it, r):
    it = list(it)
    if r == 0: return [()]
    if r > len(it): return []
    res = []
    def rec(start, cur):
        if len(cur) == r: res.append(tuple(cur)); return
        for i in range(start, len(it)):
            cur.append(it[i]); rec(i + 1, cur); cur.pop()
    rec(0, [])
    return res


class Ambient:
    def __init__(self, rs, k, dmax):
        self.pr = Proj(rs, k, dmax=dmax)
        self.rs, self.dmax = rs, dmax
        self.cols = {}
        for d in range(dmax + 1):
            for n_, c in enumerate(self.pr.basis[d]):
                self.cols.setdefault((d, c[-1]), []).append(n_)
        self.dim = {key: len(v) for key, v in self.cols.items()}

    def embed(self, d, i, vec):
        out = [0] * len(self.pr.basis[d])
        for c, a in zip(self.cols[(d, i)], vec):
            out[c] = a % P
        return out

    def arrow(self, d, i, vec, j):
        """image in P[d+1][j] coordinates."""
        full = self.pr.act(d, self.embed(d, i, vec), j)
        return [full[c] for c in self.cols.get((d + 1, j), [])]

    def preimage(self, d, i, Nnext):
        """basis of {x in P[d][i] : arrow(x,j) in Nnext[(d+1,j)] for all j}."""
        m = self.dim.get((d, i), 0)
        if m == 0: return []
        rowsA = []
        for j in self.rs.I:
            key = (d + 1, j)
            if self.dim.get(key, 0) == 0: continue
            Nb = Nnext.get(key, [])
            R, piv = rref([r[:] for r in Nb], self.dim[key])
            for bidx in range(m):
                e = [1 if t == bidx else 0 for t in range(m)]
                img = self.arrow(d, i, e, j)
                for r, c in zip(R, piv):
                    if img[c] % P:
                        g = img[c]
                        img = [(img[t] - g * r[t]) % P for t in range(len(img))]
                rowsA.append((bidx, j, img))
        # build matrix: columns = basis of P[d][i], rows = residual coords
        eqs = {}
        for bidx, j, img in rowsA:
            for t, val in enumerate(img):
                eqs.setdefault((j, t), [0] * m)
                eqs[(j, t)][bidx] = val % P
        M = [row for row in eqs.values() if any(x % P for x in row)]
        if not M:
            return [[1 if t == b else 0 for t in range(m)] for b in range(m)]
        R, piv = rref(M, m)
        free = [c for c in range(m) if c not in piv]
        ker = []
        for f in free:
            v = [0] * m
            v[f] = 1
            for r, c in zip(R, piv):
                v[c] = (-r[f]) % P
            ker.append(v)
        return ker


def all_quotient_dims(rs, k, v, dmax):
    """Every achievable graded-dimension profile for a quotient with dim vec v."""
    A = Ambient(rs, k, dmax)
    target = dict(zip(rs.I, v))
    results = {}

    def rec(d, N, acc):
        if any(acc.get(i, 0) > target.get(i, 0) for i in rs.I):
            return
        if d < 0:
            if all(acc.get(i, 0) == target.get(i, 0) for i in rs.I):
                rows = {}
                for (dd, ii), cols in A.cols.items():
                    q = len(cols) - len(N.get((dd, ii), []))
                    if q: rows.setdefault(dd, {})[ii] = q
                sig = tuple(sorted((x, tuple(sorted(y.items()))) for x, y in rows.items()))
                results[sig] = rows
            return
        verts = [i for i in rs.I if A.dim.get((d, i), 0)]
        pres = {i: A.preimage(d, i, N) for i in verts}
        opts = {i: subspaces(len(pres[i])) for i in verts}
        for combo in product(*[opts[i] for i in verts]):
            N2 = dict(N); acc2 = dict(acc)
            ok = True
            for i, sub in zip(verts, combo):
                basis = pres[i]
                rows = [[sum(coef[t] * basis[t][c] for t in range(len(basis))) % P
                         for c in range(A.dim[(d, i)])] for coef in sub]
                R, _ = rref(rows, A.dim[(d, i)])
                N2[(d, i)] = R
                q = A.dim[(d, i)] - len(R)
                acc2[i] = acc2.get(i, 0) + q
                if acc2[i] > target.get(i, 0): ok = False; break
            if ok:
                rec(d - 1, N2, acc2)

    rec(dmax, {}, {})
    return results


if __name__ == "__main__":
    from rootsys import RootSystem
    from wordtrace import trace, v_from_trace
    for label, (typ, n, k), word in [
        ("Ex 5.10", ("D",4,2), (1,2,3,4,2)),
        ("Ex 5.11", ("D",5,2), (2,1,3,4,5,3,2)),
        ("Ex 5.12", ("D",5,2), (2,1,2,3,4,5,3,2)),
    ]:
        rs = RootSystem(typ, n)
        v = v_from_trace(rs, trace(rs, word, rs.fundamental(k)))
        res = all_quotient_dims(rs, k, v, dmax=7)
        print(f"{label} v={v}: {len(res)} graded-dimension profile(s)")
        for rows in res.values():
            print("   ", {d: rows[d] for d in sorted(rows)},
                  "simple" if all(c<=1 for r in rows.values() for c in r.values()) else "NOT simple")
