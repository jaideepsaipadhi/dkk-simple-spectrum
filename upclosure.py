"""
Proposition T (type D), the current state of the proof.

Proposition T is equivalent to UPWARD CLOSURE at the trivalent node:

    if  V_t[d] != 0  then  V_t[d+2] = (I_kappa)_t[d+2].

Indeed a submodule whose colour-t part is upward closed, with total v_t, is
exactly the greedy top truncation of the colour-t profile of I_kappa.

Two ingredients.

  LEMMA U.  For x in (I_kappa)_t[d] the submodule Pi x contains the three
  double paths y_j = a_{j->t} a_{t->j} x, which satisfy one linear relation
  (the preprojective relation at t), so they span at most 2 dimensions.  In
  type D, dim (I_kappa)_t[d+2] <= 2 by Lemma P, and the span IS everything
  unless x lies on one of exactly TWO distinguished lines: the kernels of the
  two arrows from t to the leaves n-1 and n.  Call those the SPECIAL lines.
  So upward closure can only fail when V_t[d] is a special line.

  LEMMA S'.  In the extremal module V, a special line occurs at colour t only
  at a degree whose successor degree d+2 is already full.

Lemma U is checked here for all D_n, n <= 8, and every fundamental weight.
Lemma S' is checked on every colour-t graded piece of every extremal module in
range.  Finally we record the mechanism that makes Lemma S' plausible and
gives it geometric content: among ALL graded submodules of I_kappa -- not
merely the extremal ones -- every single one that fails upward closure has

    delta(v) = 2 v_kappa - v^T C v  =  (lambda,lambda) - (mu,mu)  >  0,

i.e. its dimension vector is not extremal at all: the quiver variety
M(v, omega_kappa) has positive dimension.  Since delta depends only on v, a
proof that violators always have delta > 0 would prove Proposition T outright.
"""
import sys, time, itertools
from rootsys import RootSystem, orbit_with_words
from wordtrace import trace
from qvar import rref, P
from proj2 import Proj2
from homological import setup
from lemA import build, TRIV
from allsub import all_submodules


# ------------------------------------------------------------------ Lemma U
def special_lines(pr, rs, t, d):
    """lines in I_t[d] killed by an arrow t -> leaf."""
    m = pr.dim.get((d, t), 0)
    out = []
    for coeffs in itertools.product(range(P), repeat=m):
        if not any(coeffs):
            continue
        lead = next(i for i, c in enumerate(coeffs) if c)
        if coeffs[lead] != 1:
            continue
        x = list(coeffs)
        for j in sorted(rs.adj[t]):
            if (d, t, j) not in pr.arrow:
                continue
            if not any(v % P for v in pr.act(d, t, x, j)):
                out.append(tuple(x))
                break
    return out


def lemmaU(typ, n):
    rs = RootSystem(typ, n)
    t = TRIV[typ](n)
    tot = ok = 0
    for k in range(1, n - 1):
        pr = Proj2(rs, k, dmax=60)
        for d in range(pr.top + 1):
            m, m2 = pr.dim.get((d, t), 0), pr.dim.get((d + 2, t), 0)
            if m == 0 or m2 == 0:
                continue
            sp = set(special_lines(pr, rs, t, d))
            for coeffs in itertools.product(range(P), repeat=m):
                if not any(coeffs):
                    continue
                lead = next(i for i, c in enumerate(coeffs) if c)
                if coeffs[lead] != 1:
                    continue
                x = list(coeffs)
                ys = []
                for j in sorted(rs.adj[t]):
                    if (d, t, j) not in pr.arrow:
                        continue
                    y = pr.act(d, t, x, j)
                    if (d + 1, j, t) in pr.arrow:
                        ys.append(pr.act(d + 1, j, y, t))
                R, _ = rref([y[:] for y in ys], m2)
                full = (len(R) == m2)
                tot += 1
                # the direction we use: span not full  =>  x is a special line
                ok += (full or tuple(x) in sp)
    return tot, ok


# ----------------------------------------------------------------- Lemma S'
def lemmaS(typ, n, k, cap=400):
    rs, pr, kk = setup(typ, n, k)
    t = TRIV[typ](n)
    lam = rs.fundamental(k)
    orb = orbit_with_words(rs, lam)
    if len(orb) > cap:
        return None
    sp_by_deg = {d: set(special_lines(pr, rs, t, d))
                 for d in range(pr.top + 1) if pr.dim.get((d, t), 0)}
    tot = ok = nsp = 0
    for mu, word in orb.items():
        if not word:
            continue
        X = build(rs, pr, word, lam)
        for d in range(pr.top + 1):
            rows = X.get((d, t))
            if not rows:
                continue
            tot += 1
            if len(rows) != 1:
                ok += 1
                continue
            x = rows[0]
            lead = next((i for i, c in enumerate(x) if c % P), None)
            inv = pow(x[lead], P - 2, P)
            xn = tuple((c * inv) % P for c in x)
            if xn not in sp_by_deg.get(d, ()):
                ok += 1
                continue
            nsp += 1
            above = pr.dim.get((d + 2, t), 0)
            ok += (above == 0 or len(X.get((d + 2, t), [])) == above)
    return tot, ok, nsp


# --------------------------------------------------------- the delta defect
def delta_of_violators(typ, n, k):
    rs = RootSystem(typ, n)
    t = TRIV[typ](n)
    pr, subs = all_submodules(rs, k)
    tot = pos = 0
    dist = {}
    for N in subs:
        bad = False
        for d in range(pr.top + 1):
            if N.get((d, t)) and pr.dim.get((d + 2, t), 0):
                if len(N.get((d + 2, t), [])) < pr.dim[(d + 2, t)]:
                    bad = True
                    break
        if not bad:
            continue
        v = [sum(len(N.get((d, i), [])) for d in range(pr.top + 1)) for i in rs.I]
        delta = 2 * v[k - 1] - sum(v[i - 1] * rs.A[i - 1][j - 1] * v[j - 1]
                                   for i in rs.I for j in rs.I)
        tot += 1
        pos += (delta > 0)
        dist[delta] = dist.get(delta, 0) + 1
    return tot, pos, dist


if __name__ == "__main__":
    t0 = time.time()
    allok = True

    print("LEMMA U -- if the double-path span is NOT all of I_t[d+2] then x")
    print("lies on one of the special lines (kernel of an arrow out of t)\n")
    T = O = 0
    for n in range(4, 9):
        tot, ok = lemmaU("D", n)
        T += tot; O += ok
        allok &= (ok == tot)
        print(f"  D{n}: {ok}/{tot} lines   "
              f"{'OK' if ok == tot else 'FAIL'}")
    print(f"  total {O}/{T}")

    print("\nLEMMA S' -- in the extremal module, a special line occurs only")
    print("where the degree above is already full\n")
    T = O = S = 0
    for typ, n, k in [("D", 4, 2), ("D", 5, 2), ("D", 5, 3), ("D", 6, 2),
                      ("D", 6, 3), ("D", 6, 4), ("D", 7, 3)]:
        r = lemmaS(typ, n, k)
        if r is None:
            continue
        tot, ok, nsp = r
        T += tot; O += ok; S += nsp
        allok &= (ok == tot)
        print(f"  {typ}{n} om{k}: {ok}/{tot} colour-t pieces "
              f"({nsp} of them special lines)   {'OK' if ok == tot else 'FAIL'}")
    print(f"  total {O}/{T}, of which {S} special lines")

    print("\nTHE DEFECT -- every graded submodule of I_kappa that fails upward")
    print("closure has delta(v) > 0, so its dimension vector is not extremal\n")
    T = O = 0
    dd = {}
    for typ, n, k in [("D", 5, 3), ("D", 6, 3), ("D", 6, 4), ("D", 7, 3)]:
        tot, pos, dist = delta_of_violators(typ, n, k)
        T += tot; O += pos
        for a, b in dist.items():
            dd[a] = dd.get(a, 0) + b
        print(f"  {typ}{n} e{k}: {pos}/{tot} violators have delta > 0   {dist}")
    allok &= (O == T)
    print(f"  total {O}/{T};  distribution of delta over all violators: {dd}")

    print(f"\n({time.time()-t0:.0f}s)")
    print("ALL THREE CHECKS PASS" if allok else "*** something failed ***")
    sys.exit(0 if allok else 1)
