"""
P = Pi e_k built WITHOUT enumerating walks.

Recursively:  Free_{d+1}[j] = (+)_{i ~ j} P_d[i]   (one copy per neighbour)
and P_{d+1}[j] is that modulo the preprojective relations coming from degree
d-1: for each basis vector x of P_{d-1}[j],

      sum_{i ~ j} eps(j,i) * ( arrow_{j->i}(x) placed in the i-component ) = 0.

The arrow P_d[i] -> P_{d+1}[j] is inclusion into the i-component followed by
the quotient projection.  Piece dimensions stay tiny, so this is fast at any
rank -- no 3^d blow-up.
"""
from qvar import rref, P


class Proj2:
    def __init__(self, rs, k, dmax=30):
        self.rs, self.k = rs, k
        eps = lambda a, b: 1 if a < b else -1
        self.dim = {}                    # (d,i) -> dim
        self.arrow = {}                  # (d,i,j) -> matrix rows: P_d[i] -> P_{d+1}[j]
        self.dim[(0, k)] = 1
        for i in rs.I:
            self.dim.setdefault((0, i), 0)
        d = 0
        while d < dmax:
            # free part at degree d+1, vertex j: blocks indexed by neighbours i
            blocks, free_dim = {}, {}
            for j in rs.I:
                off, tot = {}, 0
                for i in sorted(rs.adj[j]):
                    off[i] = tot
                    tot += self.dim.get((d, i), 0)
                blocks[j] = off
                free_dim[j] = tot
            # relations from degree d-1
            rel = {j: [] for j in rs.I}
            if d >= 1:
                for j in rs.I:
                    for b in range(self.dim.get((d - 1, j), 0)):
                        row = [0] * free_dim[j]
                        for i in sorted(rs.adj[j]):
                            M = self.arrow.get((d - 1, j, i))
                            if not M: continue
                            img = M[b]
                            for t, val in enumerate(img):
                                row[blocks[j][i] + t] = (row[blocks[j][i] + t]
                                                         + eps(j, i) * val) % P
                        if any(x % P for x in row):
                            rel[j].append(row)
            # quotient
            for j in rs.I:
                R, piv = rref(rel[j], free_dim[j])
                freecols = [c for c in range(free_dim[j]) if c not in piv]
                self.dim[(d + 1, j)] = len(freecols)
                # reduction map free -> quotient coordinates
                def red(vec, R=R, piv=piv, fc=freecols, n=free_dim[j]):
                    v = vec[:]
                    for r, c in zip(R, piv):
                        if v[c] % P:
                            g = v[c]
                            v = [(v[t] - g * r[t]) % P for t in range(n)]
                    return [v[c] % P for c in fc]
                for i in sorted(rs.adj[j]):
                    di = self.dim.get((d, i), 0)
                    M = []
                    for b in range(di):
                        vec = [0] * free_dim[j]
                        vec[blocks[j][i] + b] = 1
                        M.append(red(vec))
                    if di: self.arrow[(d, i, j)] = M
            if all(self.dim.get((d + 1, j), 0) == 0 for j in rs.I):
                self.top = d
                break
            d += 1
        else:
            self.top = dmax
        self.dim = {kk: vv for kk, vv in self.dim.items() if vv}

    def act(self, d, i, vec, j):
        M = self.arrow.get((d, i, j))
        if not M: return []
        n = len(M[0]) if M else 0
        out = [0] * n
        for b, a in enumerate(vec):
            if a % P == 0: continue
            for t in range(n):
                out[t] = (out[t] + a * M[b][t]) % P
        return out


if __name__ == "__main__":
    from rootsys import RootSystem
    for typ, n, k in [("D",4,2),("D",5,2),("D",6,2),("E",6,2),("D",7,2),("E",7,1)]:
        rs = RootSystem(typ, n); pr = Proj2(rs, k)
        tot = {}
        for (d, i), c in pr.dim.items(): tot[i] = tot.get(i, 0) + c
        print(f"{rs.name()} e_{k}: top degree {pr.top}, "
              f"dim vector {tuple(tot.get(i,0) for i in rs.I)}, "
              f"max piece {max(pr.dim.values())}")
