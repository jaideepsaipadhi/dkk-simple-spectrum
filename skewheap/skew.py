"""Excited diagrams and the skew hook formula, per Naruse-Okada Def 3.1 / Thm 1.2.

Def 3.1(a): u in D is D-active if there is v in (P\\D)_{c(u)} with v < u, [v,u] a
d_k-interval, and [v,u] cap D cap N_{c(u)} = empty, where N_{c} denotes the elements whose
colour is adjacent to c in the Dynkin diagram.
Def 3.1(b): the elementary excitation replaces u by the bottom element v of [v,u].
Def 3.1(c): excited diagrams of F in P = all subsets reachable from F by such moves.

Everything is in PROCTOR orientation (corner = unique maximum), the same convention Def 2.6
needs, obtained by reversing the heap order.  Colours are the heap's letters.

The univariate skew formula is tested in the form
        F(P\\F; q)  =?  sum over excited diagrams D of  prod_{u not in D} 1/(1 - q^{h(u)}),
which collapses to Peterson-Proctor when F is empty.
usage: python3 skew.py [MAXLEN=8] [PMAX=10] [TYPES=A4,D4] [MAXF=2]
"""
import sys, itertools, collections
from fractions import Fraction as Fr
sys.path.insert(0, '.')
from nsl import cartan, heap_down, is_reduced, adjacent
from dcomp2 import dcomplete
from novel import dual, relabel_natural, ser
from no26 import leq_matrix, dk_intervals, no_hooks
from filterhook import hooks_of

def to_proctor(down, word):
    """reverse the heap order; return (down_P, colours), naturally labelled."""
    n = len(down)
    up = [[] for _ in range(n)]
    for u in range(n):
        for v in down[u]: up[v].append(u)
    dP = [frozenset(up[u]) for u in range(n)]          # x <= y in P iff y <= x in heap
    need = [set(dP[u]) for u in range(n)]
    order = []; placed = set()
    while len(order) < n:
        for u in range(n):
            if u in placed: continue
            if need[u] <= placed: order.append(u); placed.add(u); break
        else: raise RuntimeError('cycle')
    pos = {u: i for i, u in enumerate(order)}
    # also return `order`: order[i] is the HEAP index of P-element i
    return ([frozenset(pos[v] for v in dP[u]) for u in order],
            [word[u] for u in order], order)

def excited_diagrams(dP, col, C, F, cap=20000):
    n = len(dP); leq = leq_matrix(dP)
    dks = {}
    for (v, u, x, y) in dk_intervals(dP, leq):
        dks.setdefault(u, []).append(v)
    inter = {}
    for u in dks:
        for v in dks[u]:
            inter[(v, u)] = [w for w in range(n) if leq[v][w] and leq[w][u]]
    seen = {frozenset(F)}; stack = [frozenset(F)]
    while stack:
        D = stack.pop()
        for u in D:
            for v in dks.get(u, []):
                if v in D: continue
                if col[v] != col[u]: continue
                I = inter[(v, u)]
                if any(w in D and adjacent(C, col[w], col[u]) for w in I): continue
                D2 = frozenset((set(D) - {u}) | {v})
                if D2 not in seen:
                    seen.add(D2); stack.append(D2)
                    if len(seen) > cap: return None
    return seen

def series_prod(hooks, N):
    f = [0]*(N+1); f[0] = 1
    for h in hooks:
        for d in range(h, N+1): f[d] += f[d-h]
    return f

def sub_series(dP, S, N):
    """generating function of (induced subposet on S)-partitions"""
    S = sorted(S); pos = {u: i for i, u in enumerate(S)}
    sp = [frozenset(pos[v] for v in dP[u] if v in pos) for u in S]
    return ser(relabel_natural(sp), N)

def main():
    ML   = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    PMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 10
    TY   = (sys.argv[3] if len(sys.argv) > 3 else "A4,D4").split(",")
    MAXF = int(sys.argv[4]) if len(sys.argv) > 4 else 2
    N = 20
    stat = collections.Counter(); ex = []
    for t in TY:
        typ, rk = t[0], int(t[1:]); C = cartan(typ, rk)
        out = []
        def rec(w):
            if w and is_reduced(C, w): out.append(tuple(w))
            if len(w) == ML: return
            for i in range(rk): rec(w + [i])
        rec([])
        for w in out:
            down = heap_down(C, list(w)); n = len(down)
            if n < 3 or n > PMAX: continue
            dc = (len(dcomplete(relabel_natural(dual(down)))) == 0)
            dP, col, _ord = to_proctor(down, list(w))
            h = no_hooks(dP)
            if h is None: continue
            leq = leq_matrix(dP)
            # order filters F of P (upper sets), small ones
            ups = [frozenset(x for x in range(n) if leq[u][x]) for u in range(n)]
            filters = set()
            for r in range(1, MAXF+1):
                for S in itertools.combinations(range(n), r):
                    F = frozenset().union(*[ups[u] for u in S])
                    if len(F) <= MAXF+2: filters.add(F)
            for F in filters:
                if not F: continue
                E = excited_diagrams(dP, col, C, F)
                if E is None: stat['too many excited diagrams'] += 1; continue
                lhs = sub_series(dP, [u for u in range(n) if u not in F], N)
                rhs = [0]*(N+1)
                for D in E:
                    pr = series_prod([h[u] for u in range(n) if u not in D], N)
                    for d in range(N+1): rhs[d] += pr[d]
                agree = lhs[:N+1] == rhs[:N+1]
                stat[(dc, agree)] += 1
                if not agree and len(ex) < 3:
                    ex.append((t, w, n, dc, sorted(F), len(E), lhs[:7], rhs[:7]))
    print(f"types {TY}, words to length {ML}, |P| <= {PMAX}, filters up to {MAXF+2} elements\n")
    for k in sorted(stat, key=str):
        print(f"   {k} -> {stat[k]}")
    for a in ex:
        print(f"\n   MISMATCH {a[0]} w={a[1]} |P|={a[2]} d-complete={a[3]} F={a[4]} "
              f"#excited={a[5]}")
        print(f"      LHS (P\\F)-partitions = {a[6]}")
        print(f"      RHS excited-diagram  = {a[7]}")

if __name__ == "__main__":
    main()
