"""
GLS (Kac-Moody groups and cluster algebras, Sect. 2.4) iterated socle:

  soc_(j_1,...,j_t)(X) = X_t  where  0 = X_0 <= X_1 <= ... <= X_t <= X
  and  X_p / X_{p-1} = soc_(j_p)( X / X_{p-1} ),
  soc_(j)(Y) = sum of all submodules of Y isomorphic to S_j.

Concretely X_p = X_{p-1} + { x in P[d][j_p] : every arrow sends x into X_{p-1} }.

V_k := soc_(i_k, ..., i_1)( I_{i_k} ) is then an EXPLICIT submodule of the
injective, built from the reduced word.  This is the construction the earlier
literature (Savage-Tingley) left abstract.
"""
from qvar import rref, P
from proj2 import Proj2


def socle_step(pr, rs, X, j):
    """X + soc_j(P/X):  add every x in P[*][j] killed into X by all arrows."""
    Y = {k: [r[:] for r in v] for k, v in X.items()}
    for d in range(pr.top + 1):
        m = pr.dim.get((d, j), 0)
        if m == 0:
            continue
        # solve for x with act(x, i) in X[(d+1,i)] for every neighbour i
        eqs = []
        for i in rs.I:
            ni = pr.dim.get((d + 1, i), 0)
            if ni == 0 or (d, j, i) not in pr.arrow:
                continue
            R, piv = rref([r[:] for r in Y.get((d + 1, i), [])], ni)
            imgs = []
            for b in range(m):
                e = [1 if x == b else 0 for x in range(m)]
                img = pr.act(d, j, e, i)
                for r, c in zip(R, piv):
                    if img[c] % P:
                        g = img[c]
                        img = [(img[x] - g * r[x]) % P for x in range(ni)]
                imgs.append(img)
            for x in range(ni):
                row = [imgs[b][x] % P for b in range(m)]
                if any(y % P for y in row):
                    eqs.append(row)
        if not eqs:
            ker = [[1 if x == b else 0 for x in range(m)] for b in range(m)]
        else:
            R, piv = rref(eqs, m)
            free = [c for c in range(m) if c not in piv]
            ker = []
            for f in free:
                x = [0] * m
                x[f] = 1
                for r, c in zip(R, piv):
                    x[c] = (-r[f]) % P
                ker.append(x)
        if ker:
            rows = Y.get((d, j), []) + ker
            Y[(d, j)] = rref(rows, m)[0]
    return Y


def iterated_socle(pr, rs, seq):
    X = {}
    for j in seq:
        X = socle_step(pr, rs, X, j)
    return X


if __name__ == "__main__":
    import json
    from rootsys import RootSystem
    from wordtrace import trace, v_from_trace
    print("Does soc_(i_k,...,i_1)(Pi e_kappa) reproduce the extremal submodule?\n")
    for typ, n, k, word in [("D",4,2,(1,2,3,4,2)),
                            ("D",5,2,(2,1,3,4,5,3,2)),
                            ("D",5,2,(2,1,2,3,4,5,3,2))]:
        rs = RootSystem(typ, n); pr = Proj2(rs, k)
        v = v_from_trace(rs, trace(rs, word, rs.fundamental(k)))
        for name, seq in [("word as given", list(word)),
                          ("word reversed", list(reversed(word)))]:
            X = iterated_socle(pr, rs, seq)
            q = tuple(sum(pr.dim.get((d,i),0) - len(X.get((d,i),[]))
                          for d in range(pr.top+1)) for i in rs.I)
            print(f"  {rs.name()} w={''.join('s'+str(i) for i in word)}  {name:<14} "
                  f"quotient dimvec {q}   target v={tuple(int(x) for x in v)}  "
                  f"{'MATCH' if q==tuple(int(x) for x in v) else ''}")
        print()
