"""
Identify the "special" colour-t vectors of I_kappa in type D: those x in
I_t[d] for which the three double paths a_{j->t}a_{t->j}x span only a line
instead of all of I_t[d+2].

For each such x we record which of the three arrows out of t kill x, i.e. the
support of x under  x |-> (a_{t->j} x)_{j ~ t}.
"""
import sys, itertools
from rootsys import RootSystem
from proj2 import Proj2
from qvar import rref, P
from homological import setup
from lemA import TRIV


def run(typ, n, k):
    rs, pr, kk = setup(typ, n, k)
    t = TRIV[typ](n)
    print(f"=== {typ}{n} om{k}   t = {t}, neighbours {sorted(rs.adj[t])}")
    for d in range(pr.top + 1):
        m = pr.dim.get((d, t), 0)
        m2 = pr.dim.get((d + 2, t), 0)
        if m == 0 or m2 == 0:
            continue
        lines = []
        for coeffs in itertools.product(range(P), repeat=m):
            if not any(coeffs) or coeffs[next(i for i, c in enumerate(coeffs) if c)] != 1:
                continue                      # one representative per line
            x = list(coeffs)
            ys, killed, zero_out = [], [], []
            for j in sorted(rs.adj[t]):
                if (d, t, j) not in pr.arrow:
                    zero_out.append(j); continue
                y = pr.act(d, t, x, j)
                if not any(v % P for v in y):
                    killed.append(j)
                if (d + 1, j, t) in pr.arrow:
                    ys.append(pr.act(d + 1, j, y, t))
            R, _ = rref([y[:] for y in ys], m2)
            if len(R) < m2:
                lines.append((tuple(x), len(R), tuple(killed), tuple(zero_out)))
        if lines:
            print(f"  degree {d}: dim I_t[d]={m}, dim I_t[d+2]={m2}; "
                  f"{len(lines)} special line(s)")
            for x, rk, killed, zo in lines:
                print(f"      x={x}  rank={rk}  arrows t->j killing x: {killed}"
                      f"  (no arrow to {zo})")


if __name__ == "__main__":
    for typ, n, k in [("D", 5, 3), ("D", 6, 3), ("D", 6, 4), ("D", 7, 5)]:
        run(typ, n, k)
        print()
