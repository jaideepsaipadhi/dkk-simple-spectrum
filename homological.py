"""
THE HOMOLOGICAL IDENTITY behind Lemma S.

Let I = I_kappa be the injective envelope of S_kappa over the preprojective
algebra Pi of a Dynkin quiver (Pi is self-injective, so I = Pi e_{nu^{-1}kappa}).
Let X <= I be a submodule with dim X = v', and put mu' = lambda - v' (root
coordinates), c = <alpha_i^vee, mu'>.

From  0 -> X -> I -> Q -> 0  and Ext^1(S_i, I) = 0 (I injective):

    dim soc_i(Q) = dim Hom(S_i,I) - dim Hom(S_i,X) + dim Ext^1(S_i,X).

Crawley-Boevey's formula for the preprojective algebra,

    dim Ext^1(M,N) = dim Hom(M,N) + dim Hom(N,M) - (dim M, dim N),

with M = S_i, N = X, and  (alpha_i, v') = 2v'_i - sum_{j~i} v'_j,  gives

    ***  dim soc_i(I/X) = c + dim top_i(X)  ***

where top_i(X) = dim Hom(X, S_i) is the multiplicity of S_i in the head of X
and c = delta_{i,kappa} + sum_{j~i} v'_j - 2 v'_i.

So a socle step of colour i enlarges X by exactly c PRECISELY WHEN X has no
S_i in its head.  This script tests all three quantities at every step.
"""
import sys, time
from rootsys import RootSystem, orbit_with_words
from wordtrace import trace, v_from_trace
from proj2 import Proj2
from gls import socle_step
from qvar import rref, P


def socv(rs, j):
    pr = Proj2(rs, j)
    return max((i for (d, i), m in pr.dim.items() if d == pr.top and m), default=None)


def setup(typ, n, k, _c={}):
    key = (typ, n, k)
    if key not in _c:
        rs = RootSystem(typ, n)
        inv = {socv(rs, j): j for j in rs.I}
        _c[key] = (rs, Proj2(rs, inv[k]), inv[k])
    return _c[key]


def dimvec(rs, X, top):
    return [sum(len(X.get((d, i), [])) for d in range(top + 1)) for i in rs.I]


def top_i(pr, rs, X, i):
    """dim Hom(X, S_i) = dim X_i - dim (sum of images of arrows j->i inside X).

    The head of X at vertex i is X_i modulo the span of arrow_{j->i}(X_j) for
    all neighbours j.  Computed degree by degree: arrows raise degree by 1, so
    the relevant image in X[(d,i)] comes from X[(d-1,j)].
    """
    tot = 0
    for d in range(pr.top + 1):
        m = pr.dim.get((d, i), 0)
        if m == 0 or not X.get((d, i)):
            continue
        img = []
        for j in rs.adj[i]:
            if (d - 1, j, i) not in pr.arrow:
                continue
            for row in X.get((d - 1, j), []):
                img.append(pr.act(d - 1, j, row, i))
        # head_i at degree d = dim X[(d,i)] - dim (X[(d,i)] cap span(img))
        # since img <= X[(d,i)] automatically (X is a submodule):
        R, _ = rref([r[:] for r in img], m)
        tot += len(X[(d, i)]) - len(R)
    return tot


def run(typ, n, k, cap=400):
    rs, pr, kk = setup(typ, n, k)
    lam = rs.fundamental(k)
    orb = orbit_with_words(rs, lam)
    if len(orb) > cap:
        return None
    steps = ok_head = ok_growth = 0
    worst = []
    for mu, word in orb.items():
        if not word:
            continue
        tr = trace(rs, word, lam)
        X = {}
        for _, i, c, _ in tr:
            v0 = dimvec(rs, X, pr.top)
            h = top_i(pr, rs, X, i)
            # predicted c from the dimension vector alone
            cpred = (1 if i == k else 0) + sum(v0[j - 1] for j in rs.adj[i]) - 2 * v0[i - 1]
            for _ in range(c):
                X = socle_step(pr, rs, X, i)
            v1 = dimvec(rs, X, pr.top)
            g = v1[i - 1] - v0[i - 1]
            steps += 1
            ok_head += (h == 0)
            ok_growth += (g == c == cpred)
            if h != 0 and len(worst) < 5:
                worst.append((typ, n, k, i, c, h))
    return steps, ok_head, ok_growth, worst


CASES = [("A", 4, 2), ("A", 5, 3), ("D", 4, 2), ("D", 5, 2), ("D", 5, 3),
         ("D", 6, 3), ("E", 6, 2), ("E", 6, 3), ("E", 6, 5), ("E", 7, 1)]

if __name__ == "__main__":
    t0 = time.time()
    S = H = G = 0
    bad = []
    for typ, n, k in CASES:
        r = run(typ, n, k)
        if r is None:
            continue
        s, h, g, w = r
        S += s; H += h; G += g; bad += w
        print(f"  {typ}{n} om{k}: head_i(X)=0 in {h}/{s} steps; "
              f"growth = c = predicted in {g}/{s}")
    print(f"\ntotals: {H}/{S} steps with no S_i in the head, "
          f"{G}/{S} with exact predicted growth  ({time.time()-t0:.0f}s)")
    if bad:
        print("counterexamples to head_i(X)=0:", bad)
    print("IDENTITY CONFIRMED" if H == S == G else "*** see above ***")
    sys.exit(0 if H == S == G else 1)
