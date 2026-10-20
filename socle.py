"""soc_sigma(I_ka) computed directly, as a linear recursion along sigma.

The socle step of colour j (typeD.tex, "Iterated socles") is
    X  |->  { x in M : a_{j->i} x in X  for every i ~ j }  (enlarging X only at vertex j),
so each step is a kernel computation, not a search.  sigma is the c-weighted sequence of a
reduced word: colour i_p repeated c_p times, in trace order.

Validated against dfs3.solve wherever that terminates; see __main__.
"""
from qvar import rref, P
from proj2 import Proj2
from wordtrace import trace

_CACHE = {}
def proj(rs, k):
    key = (rs.typ, rs.n, k)
    if key not in _CACHE: _CACHE[key] = Proj2(rs, k)
    return _CACHE[key]

def _socle_at(pr, rs, d, j, X):
    """basis (rref rows) of { x in P_d[j] : a_{j->i} x in X_{d+1,i} for all i ~ j }"""
    m = pr.dim.get((d, j), 0)
    if m == 0: return []
    eqs = []
    for i in rs.I:
        ni = pr.dim.get((d + 1, i), 0)
        if ni == 0 or (d, j, i) not in pr.arrow: continue
        R, piv = rref([r[:] for r in X.get((d + 1, i), [])], ni)
        imgs = []
        for b in range(m):
            e = [1 if t == b else 0 for t in range(m)]
            img = pr.act(d, j, e, i)
            for r, c in zip(R, piv):
                if img[c] % P:
                    g = img[c]
                    img = [(img[t] - g * r[t]) % P for t in range(ni)]
            imgs.append(img)
        for t in range(ni):
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

def soc_step(pr, rs, X, j):
    """one socle step of colour j; returns the enlarged X"""
    Y = dict(X)
    for d in range(pr.top + 1):
        if pr.dim.get((d, j), 0) == 0: continue
        B = _socle_at(pr, rs, d, j, X)
        if B:
            R, _ = rref([r[:] for r in B], pr.dim[(d, j)])
            Y[(d, j)] = R
        elif (d, j) in Y:
            pass
    return Y

def dimvec(pr, rs, X):
    v = [0] * rs.n
    for (d, i), B in X.items(): v[i - 1] += len(B)
    return tuple(v)

def rows_of(pr, rs, X):
    out = {}
    for (d, i), B in X.items():
        if B: out.setdefault(d, {})[i] = len(B)
    return out

_NAK = {}
def nakayama(rs):
    """nu with soc(Pi e_{nu(j)}) = S_j, measured rather than assumed."""
    key = (rs.typ, rs.n)
    if key in _NAK: return _NAK[key]
    nu = {}
    for j in rs.I:
        pr = proj(rs, j)
        for d in range(pr.top + 1):
            for i in rs.I:
                if _socle_at(pr, rs, d, i, {}):
                    nu[i] = j
    assert sorted(nu) == list(rs.I), nu
    _NAK[key] = nu
    return nu

def soc_sigma(rs, k, word, sigma=None):
    """soc_sigma(I_k), with I_k = Pi e_{nu(k)} the injective with socle S_k.
    sigma defaults to the c-weighted sequence of `word` acting on omega_k."""
    pr = proj(rs, nakayama(rs)[k])
    if sigma is None:
        tr = trace(rs, word, rs.fundamental(k))
        sigma = [i for _, i, c, _ in tr for _ in range(c)]
    X = {}
    for j in sigma:
        X = soc_step(pr, rs, X, j)
    return X, rows_of(pr, rs, X), dimvec(pr, rs, X)

if __name__ == "__main__":
    import time
    from rootsys import RootSystem
    from wordtrace import v_from_trace
    from dfs3 import solve
    cases = [("Ex 5.10", ("D", 4, 2), (1, 2, 3, 4, 2)),
             ("Ex 5.11", ("D", 5, 2), (2, 1, 3, 4, 5, 3, 2)),
             ("Ex 5.12", ("D", 5, 2), (2, 1, 2, 3, 4, 5, 3, 2)),
             ("ex:fcbad", ("D", 4, 2), (2, 4, 3, 1, 2)),
             ("E6 c=3",  ("E", 6, 4), (2, 4, 5, 6, 3, 4, 5, 1, 3, 4))]
    print("validating soc_sigma against dfs3.solve\n")
    for label, (typ, n, k), word in cases:
        rs = RootSystem(typ, n)
        v = v_from_trace(rs, trace(rs, word, rs.fundamental(k)))
        t = time.time(); X, rows, dv = soc_sigma(rs, k, word); ts = time.time() - t
        t = time.time(); r = solve(rs, k, v); td = time.time() - t
        drows = list(r.values())[0] if r else None
        same_dv = (dv == tuple(v))
        same_rows = (rows == drows)
        M = max(rows) if rows else 0
        rev = drows and {M - d: x for d, x in rows.items()} == drows
        print(f"{label}: v={tuple(v)}  soc dimvec={dv} {'OK' if same_dv else 'MISMATCH'}"
              f"  ({ts:.2f}s vs dfs3 {td:.2f}s)")
        print(f"    soc rows  = {dict(sorted(rows.items()))}")
        print(f"    dfs3 rows = {dict(sorted(drows.items())) if drows else None}")
        print(f"    rows equal={same_rows}  rows equal after degree reversal={bool(rev)}")
